// =============================================================================
// =============================================================================
// 01: generate and export sentinel mosaic
// =============================================================================
// =============================================================================


// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX   = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER   = 'Retamap/Retamap_GEE_Exports';
var EXPORT_TO_DRIVE = false;  // toggle: true → also export to Google Drive

var years         = ee.List.sequence(2017, 2026);
var years_toMerge = ee.List.sequence(2017, 2025);
var exportYears   = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];

var roi = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');

// =============================================================================
// 1. CLOUD SCORE+
// =============================================================================
var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

var maskLowQA = function(image) {
  return image.updateMask(image.select('cs').gte(0.6));
};

// =============================================================================
// 2. PREPROCESSING FUNCTIONS
// =============================================================================
var calculateIndices = function(image) {
  var ndyi = image.normalizedDifference(['B3', 'B2']).rename('NDYI');
  return image.addBands([ndyi]);
};

var waterMasking = function(image) {
  var waterMask = image.normalizedDifference(['B3', 'B8']).lt(0);
  return image.updateMask(waterMask);
};

var addDateBands = function(image) {
  var y = ee.Image.constant(image.date().get('year')).uint8().rename('year');
  return image.addBands([y]);
};

// =============================================================================
// 3. SENTINEL-2 COLLECTION
// =============================================================================
var elevationMask = ee.Image('USGS/SRTMGL1_003').lte(1100);

var sentinel2 = ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
  .filterBounds(roi)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 40))
  .map(function(image) {
    return image.linkCollection(csPlus, ['cs']);
  })
  .map(maskLowQA)
  .map(function(image) { return image.select(['B2', 'B3', 'B4', 'B8']); })
  .map(calculateIndices)
  .map(waterMasking)
  .map(addDateBands)
  .map(function(image) { return image.updateMask(elevationMask); });

// =============================================================================
// 4. TEMPORAL FEATURE EXTRACTION
// =============================================================================
var extractMaxNDYI_NovDec = function(year) {
  return sentinel2
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(11, 12, 'month'))
    .qualityMosaic('NDYI')
    .set('year', year);
};

var calcMedianBands_Feb = function(year) {
  return sentinel2
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(2, 2, 'month'))
    .median()
    .set('year', year);
};

var maxNDYI_NovDec_Collection  = ee.ImageCollection.fromImages(years.map(extractMaxNDYI_NovDec));
var medianBands_Feb_Collection = ee.ImageCollection.fromImages(years.map(calcMedianBands_Feb));

var mergeFebWithNovDec = function(year) {
  var novDecImage = maxNDYI_NovDec_Collection
    .filter(ee.Filter.eq('year', year)).first();
  var febImage = medianBands_Feb_Collection
    .filter(ee.Filter.eq('year', ee.Number(year).add(1))).first();
  var febBands = febImage
    .select(['B2', 'B3', 'B4', 'B8', 'NDYI'])
    .rename(['B2_feb', 'B3_feb', 'B4_feb', 'B8_feb', 'NDYI_feb']);
  return novDecImage.addBands(febBands).set('year', year);
};

var mergedCollection = ee.ImageCollection(years_toMerge.map(mergeFebWithNovDec));

// =============================================================================
// 5. EXPORT mergedCollection TO ASSET
// Final band schema: B2, B3, B4, B8, NDYI, year, B2_feb, B3_feb, B4_feb, B8_feb, NDYI_feb
// Asset name prefix: 01_
// =============================================================================
exportYears.forEach(function(year) {
  var img = mergedCollection.filter(ee.Filter.eq('year', year)).first();

  Export.image.toAsset({
    image      : img,
    description: 'Export_01_MergedBands_' + year,
    assetId    : ASSET_PREFIX + '01_MergedBands_' + year,
    region     : roi,
    scale      : 10,
    maxPixels  : 1e13
  });

  if (EXPORT_TO_DRIVE) {
    Export.image.toDrive({
      image         : img,
      description   : 'Drive_01_MergedBands_' + year,
      folder        : DRIVE_FOLDER,
      fileNamePrefix: '01_MergedBands_' + year,
      region        : roi,
      scale         : 10,
      maxPixels     : 1e13,
      fileFormat    : 'GeoTIFF'
    });
  }
});
