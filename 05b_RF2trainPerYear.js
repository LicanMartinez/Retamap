// =============================================================================
// 05b: RF2 — Per-year train + predict (year-specific classifiers)
// -----------------------------------------------------------------------------
// Prerequisite : assets 04_RF2_samples_YYYY{SUFFIX}_train  (from 04_RF2sampleExport.js)
//              + assets 01_MergedBands_YYYY                 (from 01_sentinelMosaic.js)
// Produces     : assets 05b_RF2_perYear_raw_YYYY{SUFFIX}   (uint8, no patch filter)
// -----------------------------------------------------------------------------
// Unlike the main pipeline (05 + 06) which trains ONE classifier on all years
// merged, this script trains a SEPARATE RF2 per year so each model learns the
// spectral signature of its specific year. Raw predictions only — no gap-fill
// or patch filter applied.
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER    = 'Retamap/Retamap_GEE_Exports';
var EXPORT_TO_DRIVE = false;

var exportYears = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var roi         = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');
var BANDS       = ['B2', 'B3', 'B4', 'B8', 'NDYI', 'B2_feb', 'B3_feb', 'B4_feb', 'B8_feb', 'NDYI_feb'];
var RUN_SUFFIX  = '_trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k';

// =============================================================================
// 1. PER-YEAR TRAIN + PREDICT
// -----------------------------------------------------------------------------
// For each year: load that year's _train samples, train RF2(500), classify
// that year's 01_MergedBands image. Each year gets its own independent model.
// =============================================================================
var predictions = {};

exportYears.forEach(function(year) {
  var samples = ee.FeatureCollection(
    ASSET_PREFIX + '04_RF2_samples_' + year + RUN_SUFFIX + '_train'
  );

  print(year + ' — class counts:', samples.aggregate_histogram('stable_label'));

  var rf2 = ee.Classifier.smileRandomForest({numberOfTrees: 500}).train({
    features       : samples,
    classProperty  : 'stable_label',
    inputProperties: BANDS
  });

  print(year + ' — variable importance:', ee.Dictionary(rf2.explain().get('importance')));

  var img = ee.Image(ASSET_PREFIX + '01_MergedBands_' + year).select(BANDS);
  predictions[year] = img.classify(rf2).rename('classification');
});

// =============================================================================
// 2. VISUALIZATION
// =============================================================================
Map.centerObject(roi, 10);

exportYears.forEach(function(year) {
  Map.addLayer(
    predictions[year].clip(roi),
    {min: 0, max: 1, palette: ['#d3d3d3', '#2d9c00']},
    'RF2 perYear raw ' + year, false
  );
});

// =============================================================================
// 3. EXPORT RAW PREDICTIONS (05b_RF2_perYear_raw_YYYY{SUFFIX})
// =============================================================================
exportYears.forEach(function(year) {
  var img = predictions[year].select('classification').clip(roi).uint8();

  Export.image.toAsset({
    image      : img,
    description: 'Export_05b_RF2_perYear_raw_' + year + RUN_SUFFIX,
    assetId    : ASSET_PREFIX + '05b_RF2_perYear_raw_' + year + RUN_SUFFIX,
    region     : roi,
    scale      : 10,
    maxPixels  : 1e13
  });

  if (EXPORT_TO_DRIVE) {
    Export.image.toDrive({
      image         : img,
      description   : 'Drive_05b_RF2_perYear_raw_' + year + RUN_SUFFIX,
      folder        : DRIVE_FOLDER,
      fileNamePrefix: '05b_RF2_perYear_raw_' + year + RUN_SUFFIX,
      region        : roi,
      scale         : 10,
      maxPixels     : 1e13,
      fileFormat    : 'GeoTIFF'
    });
  }
});
