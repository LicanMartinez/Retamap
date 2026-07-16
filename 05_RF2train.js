// =============================================================================
// 05: RF2 — Train classifier from per-year _train sample assets → export classifier
// -----------------------------------------------------------------------------
// Prerequisite : assets 04_RF2_samples_YYYY{SUFFIX}_train  (from 04_RF2sampleExport.js)
// Produces     : asset 05_RF2_classifier{SUFFIX}
// -----------------------------------------------------------------------------
// Exporting the trained classifier as an asset means the downstream prediction
// script (06_RF2predict.js) loads it with ee.Classifier.load() instead of
// re-training. Each export task in 06 pays only the cost of classification.
// Accuracy is assessed separately by the independent validation family (09v).
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX = 'projects/ee-licanemartinez/assets/Retamap/';
var exportYears  = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var BANDS        = ['B2', 'B3', 'B4', 'B8', 'NDYI', 'B2_feb', 'B3_feb', 'B4_feb', 'B8_feb', 'NDYI_feb'];
var RUN_SUFFIX   = '_trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k';

// =============================================================================
// 1. LOAD AND MERGE PER-YEAR TRAIN SAMPLE ASSETS
// -----------------------------------------------------------------------------
// 04_RF2sampleExport.js exports one _train asset per year.
// Flatten merges them into a single FeatureCollection for training.
// =============================================================================
var samples_rf2 = ee.FeatureCollection(
  exportYears.map(function(y) {
    return ee.FeatureCollection(ASSET_PREFIX + '04_RF2_samples_' + y + RUN_SUFFIX + '_train');
  })
).flatten();

print('Class counts (stable_label):', samples_rf2.aggregate_histogram('stable_label'));

// =============================================================================
// 2. TRAIN RF2 (train samples only)
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
// 3. EXPORT CLASSIFIER AS ASSET (05_RF2_classifier{SUFFIX})
// -----------------------------------------------------------------------------
// ee.Classifier can only be exported to Asset (no Drive equivalent).
// =============================================================================
Export.classifier.toAsset({
  classifier : rf2,
  description: 'Export_05_RF2_classifier' + RUN_SUFFIX,
  assetId    : ASSET_PREFIX + '05_RF2_classifier' + RUN_SUFFIX
});
