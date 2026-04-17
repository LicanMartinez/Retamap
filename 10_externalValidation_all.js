// =============================================================================
// 10_all: Validation — external validation 2023 & Cross-run comparison (all runs)
// -----------------------------------------------------------------------------
// Prerequisite : asset 02_RF1_raw_prediction_2023              (from 02_RF1fit.js)
//              + assets 06_RF2_raw_prediction_2023{SUFFIX}     (from 06_RF2predict.js)
//              + assets 07_RF2_gapFill_2023{SUFFIX}            (from 07_RF2gapFill.js)
//              + assets 08_RF2_prediction_2023{SUFFIX}         (from 08_RF2patchFilter.js)
//              + asset 4-Validation_points_complete_2023
// Produces     : console output with full metrics for all models and runs
//              + Export.table.toDrive (CSV) for downstream analysis in R
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER    = 'Retamap/Retamap_GEE_Exports';
var EXPORT_TO_DRIVE = false;  // toggle: true → also export to Google Drive

var VAL_YEAR = 2023;

// Toggle to include or exclude the final RF2 prediction (after patch filter)
var INCLUDE_RF2_FINAL   = false; // Set to true to include 08_RF2_prediction
var INCLUDE_RF2_GAPFILL = false; // Set to true to include 07_RF2_gapFill

// List all runs to compare (ensure the primary run is included here)
var COMPARE_SUFFIXES = [
  '_s1.10k_qs1.0_s0.20k_qs0.0',
];

// =============================================================================
// 1. LOAD ASSETS
// =============================================================================

// External validation points
var valPointsRaw = ee.FeatureCollection(ASSET_PREFIX + '4-Validation_points_complete_2023');

// Convert TIPO → numeric label (retama=1, ctrl=0) for errorMatrix
var valPoints = valPointsRaw.map(function(f) {
  var label = ee.Algorithms.If(ee.String(f.get('TIPO')).equals('retama'), 1, 0);
  return f.set('label', label);
});

// =============================================================================
// 2. HELPER: build an ee.Feature with all metrics from an errorMatrix
// =============================================================================
var metricsFeature = function(label, cm, count) {
  var pa1 = ee.Array(cm.producersAccuracy()).get([1, 0]);  // recall    class-1
  var ca1 = ee.Array(cm.consumersAccuracy()).get([0, 1]);  // precision class-1
  var f1  = pa1.multiply(ca1).multiply(2).divide(pa1.add(ca1));
  return ee.Feature(null, {
    'Model'         : label,
    'Matched_Points': count,
    'OA'            : cm.accuracy(),
    'Kappa'         : cm.kappa(),
    'Recall'        : pa1,
    'Precision'     : ca1,
    'F1'            : f1
  });
};

// =============================================================================
// 3. UNIFIED VALIDATION & COMPARISON
// =============================================================================
print('═════════════════════════════════════');
print('Validation & Cross-run Comparison ' + VAL_YEAR);

var unifiedFeatures = [];

// --- 3a. RF1 (Run-independent) ---
var rf1_pred = ee.Image(ASSET_PREFIX + '02_RF1_raw_prediction_' + VAL_YEAR).select('pred');

var rf1_sampled = rf1_pred.sampleRegions({
  collection: valPoints,
  properties: ['label'],
  scale     : 10,
  tileScale : 4
});
print('RF1 — matched points:', rf1_sampled.size());

var cm_rf1 = rf1_sampled.errorMatrix({
  actual   : 'label',
  predicted: 'pred',
  order    : [0, 1]
});
unifiedFeatures.push(metricsFeature('RF1', cm_rf1, rf1_sampled.size()));

if (EXPORT_TO_DRIVE) {
  Export.table.toDrive({
    collection    : rf1_sampled,
    description   : 'Validation_RF1_external_' + VAL_YEAR,
    folder        : DRIVE_FOLDER,
    fileNamePrefix: 'Validation_RF1_external_' + VAL_YEAR,
    fileFormat    : 'CSV'
  });
}

// --- 3b. RF2 (Loop over all comparison suffixes) ---
COMPARE_SUFFIXES.forEach(function(suffix) {

  // RF2 raw
  var rf2_raw = ee.Image(ASSET_PREFIX + '06_RF2_raw_prediction_' + VAL_YEAR + suffix).select('classification');
  var rf2_raw_sampled = rf2_raw.sampleRegions({
    collection: valPoints,
    properties: ['label'],
    scale     : 10,
    tileScale : 4
  });

  var cm_rf2_raw = rf2_raw_sampled.errorMatrix({
    actual   : 'label',
    predicted: 'classification',
    order    : [0, 1]
  });
  unifiedFeatures.push(metricsFeature('RF2 raw (' + suffix + ')', cm_rf2_raw, rf2_raw_sampled.size()));

  if (EXPORT_TO_DRIVE) {
    Export.table.toDrive({
      collection    : rf2_raw_sampled,
      description   : 'Validation_RF2raw_external_' + VAL_YEAR + suffix,
      folder        : DRIVE_FOLDER,
      fileNamePrefix: 'Validation_RF2raw_external_' + VAL_YEAR + suffix,
      fileFormat    : 'CSV'
    });
  }

  // RF2 gap-fill
  if (INCLUDE_RF2_GAPFILL) {
    var rf2_gf = ee.Image(ASSET_PREFIX + '07_RF2_gapFill_' + VAL_YEAR + suffix).select('classification');
    var rf2_gf_sampled = rf2_gf.sampleRegions({
      collection: valPoints,
      properties: ['label'],
      scale     : 10,
      tileScale : 4
    });

    var cm_rf2_gf = rf2_gf_sampled.errorMatrix({
      actual   : 'label',
      predicted: 'classification',
      order    : [0, 1]
    });
    unifiedFeatures.push(metricsFeature('RF2 gap-fill (' + suffix + ')', cm_rf2_gf, rf2_gf_sampled.size()));

    if (EXPORT_TO_DRIVE) {
      Export.table.toDrive({
        collection    : rf2_gf_sampled,
        description   : 'Validation_RF2gapFill_external_' + VAL_YEAR + suffix,
        folder        : DRIVE_FOLDER,
        fileNamePrefix: 'Validation_RF2gapFill_external_' + VAL_YEAR + suffix,
        fileFormat    : 'CSV'
      });
    }
  }

  // RF2 final
  if (INCLUDE_RF2_FINAL) {
    var rf2_fin = ee.Image(ASSET_PREFIX + '08_RF2_prediction_' + VAL_YEAR + suffix).select('classification');
    var rf2_fin_sampled = rf2_fin.sampleRegions({
      collection: valPoints,
      properties: ['label'],
      scale     : 10,
      tileScale : 4
    });

    var cm_rf2_fin = rf2_fin_sampled.errorMatrix({
      actual   : 'label',
      predicted: 'classification',
      order    : [0, 1]
    });
    unifiedFeatures.push(metricsFeature('RF2 final (' + suffix + ')', cm_rf2_fin, rf2_fin_sampled.size()));

    if (EXPORT_TO_DRIVE) {
      Export.table.toDrive({
        collection    : rf2_fin_sampled,
        description   : 'Validation_RF2final_external_' + VAL_YEAR + suffix,
        folder        : DRIVE_FOLDER,
        fileNamePrefix: 'Validation_RF2final_external_' + VAL_YEAR + suffix,
        fileFormat    : 'CSV'
      });
    }
  }
});

// --- 3c. Unified Summary Table ---
print(
  ui.Chart.feature.byFeature(
    ee.FeatureCollection(unifiedFeatures), 'Model', ['Matched_Points', 'OA', 'Kappa', 'Recall', 'Precision', 'F1']
  )
    .setChartType('Table')
    .setOptions({title: 'Unified Validation & Comparison (' + VAL_YEAR + ')'})
);
