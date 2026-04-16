// =============================================================================
// 05: RF2 — Train classifier from per-year sample assets → export classifier
// -----------------------------------------------------------------------------
// Prerequisite : assets RF2_samples_YYYY  (from 04_RF2sampleExport.js)
// Produces     : asset RF2_classifier
// -----------------------------------------------------------------------------
// Exporting the trained classifier as an asset means the downstream prediction
// script (06_RF2predict.js) loads it with ee.Classifier.load() instead of
// re-training. Each export task in 06 pays only the cost of classification,
// not re-training.
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX = 'projects/ee-licanemartinez/assets/Retamap/';
var exportYears  = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var BANDS        = ['B2', 'B3', 'B4', 'B8', 'NDYI', 'B2_feb', 'B3_feb', 'B4_feb', 'B8_feb', 'NDYI_feb'];
var RUN_SUFFIX   = '_QS_n10k';

// =============================================================================
// 1. LOAD AND MERGE PER-YEAR SAMPLE ASSETS
// -----------------------------------------------------------------------------
// 04_RF2sampleExport.js exports one asset per year.
// Flatten merges them into a single FeatureCollection for training.
// =============================================================================
var samples_rf2 = ee.FeatureCollection(
  exportYears.map(function(y) {
    return ee.FeatureCollection(ASSET_PREFIX + 'RF2_samples_' + y + RUN_SUFFIX);
  })
).flatten();

// =============================================================================
// 2. TRAIN RF2 (all samples — validation is handled externally in 08)
// =============================================================================
print('Train samples:', samples_rf2.size());

var rf2 = ee.Classifier.smileRandomForest({numberOfTrees: 500}).train({
  features       : samples_rf2,
  classProperty  : 'stable_label',
  inputProperties: BANDS
});

var dict = rf2.explain();
var variableImportance = ee.Dictionary(dict.get('importance'));
print('variableImportance:', variableImportance);

// =============================================================================
// 3. EXPORT CLASSIFIER AS ASSET
// -----------------------------------------------------------------------------
// ee.Classifier can only be exported to Asset (no Drive equivalent).
// =============================================================================
Export.classifier.toAsset({
  classifier : rf2,
  description: 'Export_RF2_classifier' + RUN_SUFFIX,
  assetId    : ASSET_PREFIX + 'RF2_classifier' + RUN_SUFFIX
});
