// =============================================================================
// 08: Validation — external validation 2023 & Cross-run comparison
// -----------------------------------------------------------------------------
// Prerequisite : asset RF1_prediction_2023         (from 02_RF1fit.js)
//              + asset RF2_raw_prediction_2023     (from 06_RF2predict.js)
//              + asset RF2_prediction_2023         (from 07_RF2patchFilter.js)
//              + asset 4-Validation_points_complete_2023
// Produces     : console output with full metrics for all models and runs
//              + Export.table.toDrive (CSV) for downstream analysis in R
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER = 'Retamap/Retamap_GEE_Exports';
var VAL_YEAR     = 2023;

// List all runs to compare (ensure the primary run is included here)
var COMPARE_SUFFIXES = [
  '_noQS_n20k_compSamp',
  '_noQS1_n20k',
  // '_QS_n10k_compSamp', 
  // "_QS_n10k",
];

// =============================================================================
// 1. LOAD ASSETS
// =============================================================================

// --- 1a. External validation points ---
var valPointsRaw = ee.FeatureCollection(ASSET_PREFIX + '4-Validation_points_complete_2023');

// Convert TIPO → numeric label (retama=1, ctrl=0) for errorMatrix
var valPoints = valPointsRaw.map(function(f) {
  var label = ee.Algorithms.If(ee.String(f.get('TIPO')).equals('retama'), 1, 0);
  return f.set('label', label);
});

// =============================================================================
// 2. HELPER: build an ee.Feature with all metrics from an errorMatrix
// =============================================================================
var metricsFeature = function(label, cm) {
  var pa1 = ee.Array(cm.producersAccuracy()).get([1, 0]);  // recall    class-1
  var ca1 = ee.Array(cm.consumersAccuracy()).get([0, 1]);  // precision class-1
  var f1  = pa1.multiply(ca1).multiply(2).divide(pa1.add(ca1));
  return ee.Feature(null, {
    'Model'    : label,
    'OA'       : cm.accuracy(),
    'Kappa'    : cm.kappa(),
    'Recall'   : pa1,
    'Precision': ca1,
    'F1'       : f1
  });
};

// =============================================================================
// 3. UNIFIED VALIDATION & COMPARISON
// =============================================================================
print('═════════════════════════════════════');
print('Validation & Cross-run Comparison ' + VAL_YEAR);

var unifiedFeatures = [];

// --- 3a. RF1 (Run-independent) ---
var rf1_pred = ee.Image(ASSET_PREFIX + 'RF1_prediction_' + VAL_YEAR).select('pred');

var rf1_sampled = rf1_pred.sampleRegions({
  collection : valPoints,
  properties : ['label'],
  scale      : 10,
  tileScale  : 4
});
print('RF1 — matched points:', rf1_sampled.size());

var cm_rf1 = rf1_sampled.errorMatrix({
  actual   : 'label',
  predicted: 'pred',
  order    : [0, 1]
});
unifiedFeatures.push(metricsFeature('RF1', cm_rf1));

Export.table.toDrive({
  collection     : rf1_sampled,
  description    : 'Validation_RF1_external_' + VAL_YEAR,
  folder         : DRIVE_FOLDER,
  fileNamePrefix : 'Validation_RF1_external_' + VAL_YEAR,
  fileFormat     : 'CSV'
});

// --- 3b. RF2 (Loop over all comparison suffixes) ---
COMPARE_SUFFIXES.forEach(function(suffix) {
  
  // RF2 raw
  var rf2_raw = ee.Image(ASSET_PREFIX + 'RF2_raw_prediction_' + VAL_YEAR + suffix).select('classification');
  var rf2_raw_sampled = rf2_raw.sampleRegions({
    collection : valPoints,
    properties : ['label'],
    scale      : 10,
    tileScale  : 4
  });
  
  var cm_rf2_raw = rf2_raw_sampled.errorMatrix({
    actual   : 'label',
    predicted: 'classification',
    order    : [0, 1]
  });
  unifiedFeatures.push(metricsFeature('RF2 raw (' + suffix + ')', cm_rf2_raw));

  Export.table.toDrive({
    collection     : rf2_raw_sampled,
    description    : 'Validation_RF2raw_external_' + VAL_YEAR + suffix,
    folder         : DRIVE_FOLDER,
    fileNamePrefix : 'Validation_RF2raw_external_' + VAL_YEAR + suffix,
    fileFormat     : 'CSV'
  });

  // RF2 final
  var rf2_fin = ee.Image(ASSET_PREFIX + 'RF2_prediction_' + VAL_YEAR + suffix).select('classification');
  var rf2_fin_sampled = rf2_fin.sampleRegions({
    collection : valPoints,
    properties : ['label'],
    scale      : 10,
    tileScale  : 4
  });
  
  var cm_rf2_fin = rf2_fin_sampled.errorMatrix({
    actual   : 'label',
    predicted: 'classification',
    order    : [0, 1]
  });
  unifiedFeatures.push(metricsFeature('RF2 final (' + suffix + ')', cm_rf2_fin));

  Export.table.toDrive({
    collection     : rf2_fin_sampled,
    description    : 'Validation_RF2final_external_' + VAL_YEAR + suffix,
    folder         : DRIVE_FOLDER,
    fileNamePrefix : 'Validation_RF2final_external_' + VAL_YEAR + suffix,
    fileFormat     : 'CSV'
  });
});

// --- 3c. Unified Summary Table ---
print(
  ui.Chart.feature.byFeature(
    ee.FeatureCollection(unifiedFeatures), 'Model', ['OA', 'Kappa', 'Recall', 'Precision', 'F1']
  )
    .setChartType('Table')
    .setOptions({title: 'Unified Validation & Comparison (' + VAL_YEAR + ')'})
);