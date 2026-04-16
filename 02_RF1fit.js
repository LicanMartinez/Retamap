// =============================================================================
// 02: RF1 training, prediction, and export to asset
// Runs after 01_generate_sentinel_mosaic.js assets are exported
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER = 'Retamap/Retamap_GEE_Exports';

var exportYears = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];

var roi          = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');
var gt_polygons_new = ee.FeatureCollection(ASSET_PREFIX + 'gt_polys_6_ctrlReduce_moreCtrls_refineRetamas');  

// Band schema from phase 1:
// B2, B3, B4, B8, NDYI, B2_feb, B3_feb, B4_feb, B8_feb, NDYI_feb
var BANDS = ['B2', 'B3', 'B4', 'B8', 'NDYI', 'B2_feb', 'B3_feb', 'B4_feb', 'B8_feb', 'NDYI_feb'];

// =============================================================================
// 1. LOAD mergedCollection FROM ASSETS
// =============================================================================
var mergedCollection = ee.ImageCollection.fromImages(
  exportYears.map(function(year) {
    return ee.Image(ASSET_PREFIX + 'MergedBands_' + year).set('year', year);
  })
);

// =============================================================================
// 2. LABEL FUNCTION
// =============================================================================
var addType01 = function(feature) {
  var label = ee.Algorithms.If(
    ee.String(feature.get('TIPO')).equals('retama'), 1, 0
  );
  return feature.set('type_01', label);
};

// =============================================================================
// 3. STRATIFIED SAMPLING (two-step: fixed N points → per-year pixel lookup)
// =============================================================================
var N_CTRL   = 15000;  // total class-0 points across all years
var N_RETAMA = 10000;  // total class-1 points across all years

// Step 1: reduce ALL polygons (all years) to a 2-band raster: type_01 + ANIO_FL.
// stratifiedSample runs ONCE over the full polygon set → total N is controlled.
var labeledPolygons = gt_polygons_new.map(addType01);

var imgType = labeledPolygons.reduceToImage({
  properties: ['type_01'],
  reducer: ee.Reducer.first()
}).rename('type_01');

var imgAnioFL = labeledPolygons.reduceToImage({
  properties: ['ANIO_FL'],
  reducer: ee.Reducer.first()
}).rename('ANIO_FL');

var labelImage = imgType.addBands(imgAnioFL);

var samplePoints = labelImage
  .stratifiedSample({
    numPoints  : 0,
    classBand  : 'type_01',
    region     : roi,
    scale      : 10,
    classValues: [0, 1],
    classPoints: [N_CTRL, N_RETAMA],
    geometries : true,
    tileScale  : 16
  });
// Each point carries type_01 and ANIO_FL.
// Distribution across years is proportional to polygon area per year/class.
// Inspect with: print(samplePoints.aggregate_histogram('ANIO_FL'));

// Step 2: for each year's image, look up spectral values only at that year's points.
// sampleRegions is a point-in-raster lookup — cheap, runs in parallel across years.
var samples_rf1 = mergedCollection.map(function(image) {
  var year = image.get('year');

  var yearPoints = samplePoints.filter(ee.Filter.eq('ANIO_FL', year));

  return image.select(BANDS)
    .sampleRegions({
      collection : yearPoints,
      properties : ['type_01', 'ANIO_FL'],
      scale      : 10,
      geometries : false
    })
    .map(function(f) { return f.set('year', year); });

}).flatten();

print(samples_rf1.size());

// =============================================================================
// 4. TRAIN RF1
// =============================================================================
var rf1 = ee.Classifier.smileRandomForest({numberOfTrees: 500}).train({
  features: samples_rf1,
  classProperty: 'type_01',
  inputProperties: BANDS
});

// =============================================================================
// 5. PREDICT + CONNECTED PIXEL FILTER
// =============================================================================
var predictRF1 = function(image) {
  var pred = image.select(BANDS).classify(rf1).rename('pred');

  // // Connected pixel filter applied only to class-1 patches.
  // // selfMask() makes class-0 pixels nodata temporarily so connectedPixelCount
  // // counts only retama (1) pixels. connected is nodata where pred=0.
  // var connected = pred.selfMask().connectedPixelCount(50, true);  // 8-connected

  // // Small class-1 patches (< 5 pixels) are reclassified to 0 instead of
  // // becoming nodata. Class-0 pixels are preserved as-is.
  // // connected.unmask(0): for class-0 pixels (where connected=nodata), unmask
  // // gives 0, but pred.eq(1) is false there, so the where() condition is false
  // // and class-0 pixels are untouched.
  // var predClean = pred
  //   .where(
  //     pred.eq(1).and(connected.unmask(0).lt(5)),
  //     0
  //   )
  //   .rename('pred');

  // return predClean.set('year', image.get('year'));
  
  return pred.set('year', image.get('year'));
  
};

var preds_collection = mergedCollection.map(predictRF1);

// =============================================================================
// 6. EXPORT RF1 PREDICTIONS TO ASSET
// =============================================================================
exportYears.forEach(function(year) {
  var img = ee.Image(
    preds_collection.filter(ee.Filter.eq('year', year)).first()
  ).select('pred').clip(roi).uint8();

  Export.image.toAsset({
    image: img,
    description: 'Export_RF1_raw' + year,
    assetId: ASSET_PREFIX + 'RF1_raw_prediction_' + year,
    region: roi,
    scale: 10,
    maxPixels: 1e13
  });

  Export.image.toDrive({
    image         : img,
    description   : 'Drive_RF1_' + year,
    folder        : DRIVE_FOLDER,
    fileNamePrefix: 'RF1_raw_prediction_' + year,
    region        : roi,
    scale         : 10,
    maxPixels     : 1e13,
    fileFormat    : 'GeoTIFF'
  });
});