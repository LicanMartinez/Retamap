// =============================================================================
// 06: RF2 — Load classifier asset → raw per-year predictions
// -----------------------------------------------------------------------------
// Prerequisite : asset RF2_classifier         (from 05_RF2train.js)
//              + assets MergedBands_YYYY       (from 01_sentinelMosaic.js)
// Produces     : assets RF2_raw_prediction_YYYY  (uint8, no patch filter)
// -----------------------------------------------------------------------------
// Loading the classifier from an asset means each export task pays only the
// cost of classification, not re-training. Proceed to 07_RF2patchFilter.js
// once these assets are ingested.
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER = 'Retamap/Retamap_GEE_Exports';
var exportYears  = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var roi          = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');
var BANDS        = ['B2', 'B3', 'B4', 'B8', 'NDYI', 'B2_feb', 'B3_feb', 'B4_feb', 'B8_feb', 'NDYI_feb'];
var RUN_SUFFIX   = '_QS_n10k_compSamp';

// =============================================================================
// 1. LOAD ASSETS
// =============================================================================
var mergedCollection = ee.ImageCollection.fromImages(
  exportYears.map(function(year) {
    return ee.Image(ASSET_PREFIX + 'MergedBands_' + year).set('year', year);
  })
);

// Pre-trained classifier — no re-training per export task.
var rf2 = ee.Classifier.load(ASSET_PREFIX + 'RF2_classifier' + RUN_SUFFIX);

// =============================================================================
// 2. RAW PREDICTION PER YEAR (no patch filter)
// -----------------------------------------------------------------------------
// connectedPixelCount intentionally omitted — applied in 07_RF2patchFilter.js
// on already-materialised rasters (much cheaper).
// =============================================================================
var rf2_raw_collection = mergedCollection.map(function(image) {
  return image.select(BANDS)
    .classify(rf2)
    .rename('classification')
    .set('year', image.get('year'));
});

// =============================================================================
// 3. VISUALIZATION
// =============================================================================
Map.centerObject(roi, 10);

exportYears.forEach(function(year) {
  Map.addLayer(
    ee.Image(rf2_raw_collection.filter(ee.Filter.eq('year', year)).first()).clip(roi),
    {min: 0, max: 1, palette: ['#d3d3d3', '#2d9c00']},
    'RF2 raw ' + year, false
  );
});

// =============================================================================
// 4. EXPORT RAW PREDICTIONS PER YEAR
// =============================================================================
exportYears.forEach(function(year) {
  var img = ee.Image(
    rf2_raw_collection.filter(ee.Filter.eq('year', year)).first()
  ).select('classification').clip(roi).uint8();

  Export.image.toAsset({
    image      : img,
    description: 'Export_RF2raw_' + year + RUN_SUFFIX,
    assetId    : ASSET_PREFIX + 'RF2_raw_prediction_' + year + RUN_SUFFIX,
    region     : roi,
    scale      : 10,
    maxPixels  : 1e13
  });

  Export.image.toDrive({
    image          : img,
    description    : 'Drive_RF2raw_' + year + RUN_SUFFIX,
    folder         : DRIVE_FOLDER,
    fileNamePrefix : 'RF2_raw_prediction_' + year + RUN_SUFFIX,
    region         : roi,
    scale          : 10,
    maxPixels      : 1e13,
    fileFormat     : 'GeoTIFF'
  });
});
