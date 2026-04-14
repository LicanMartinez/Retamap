// =============================================================================
// 08: Validation — external validation 2023
// -----------------------------------------------------------------------------
// Prerequisite : asset RF1_prediction_2023         (from 02_RF1fit.js)
//              + asset RF2_raw_prediction_2023     (from 06_RF2predict.js)
//              + asset RF2_prediction_2023         (from 07_RF2patchFilter.js)
//              + asset 4-Validation_points_complete_2023
// Produces     : console output with full metrics for all three classifiers
//              + Export.table.toDrive (3 CSV) for downstream analysis in R
// -----------------------------------------------------------------------------
// Validation structure
// └─ B. External validation 2023 (independent points)
//       ├─ B1. RF1_prediction_2023
//       ├─ B2. RF2_raw_prediction_2023   (before patch filter)
//       └─ B3. RF2_prediction_2023       (after patch filter)
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER = 'Retamap/Retamap_GEE_Exports';
var VAL_YEAR     = 2023;
var RUN_SUFFIX   = '_noQS_n20k_compSamp';

// Section C — cross-run comparison (edit to include the runs you want to compare)
var COMPARE_SUFFIXES = [
  '_noQS_n20k_compSamp',
  // '_noQS_n1000',       // example: uncomment / add other runs here
];

// =============================================================================
// 1. LOAD ASSETS
// =============================================================================

// --- 1a. Prediction rasters (VAL_YEAR only) ---
// RF1 uses band 'pred'; RF2 scripts use band 'classification'
var rf1_pred     = ee.Image(ASSET_PREFIX + 'RF1_prediction_'     + VAL_YEAR).select('pred');
var rf2_raw_pred = ee.Image(ASSET_PREFIX + 'RF2_raw_prediction_' + VAL_YEAR + RUN_SUFFIX).select('classification');
var rf2_pred     = ee.Image(ASSET_PREFIX + 'RF2_prediction_'     + VAL_YEAR + RUN_SUFFIX).select('classification');

// --- 1b. External validation points ---
// Schema: TIPO ('retama'/'ctrl'), YEAR (2023)
var valPointsRaw = ee.FeatureCollection(ASSET_PREFIX + '4-Validation_points_complete_2023');

// Convert TIPO → numeric label (retama=1, ctrl=0) for errorMatrix
var valPoints = valPointsRaw.map(function(f) {
  var label = ee.Algorithms.If(ee.String(f.get('TIPO')).equals('retama'), 1, 0);
  return f.set('label', label);
});

print('External validation points (total):', valPoints.size());
print('  retama:', valPoints.filter(ee.Filter.eq('label', 1)).size());
print('  ctrl:  ', valPoints.filter(ee.Filter.eq('label', 0)).size());

// =============================================================================
// 2. HELPER: compute and print full metrics from an errorMatrix
// =============================================================================
var printMetrics = function(label, cm) {
  // producersAccuracy: ee.Array shape [[pa0], [pa1]]
  // consumersAccuracy: ee.Array shape [[ca0, ca1]]
  var pa1 = ee.Array(cm.producersAccuracy()).get([1, 0]);  // recall    class-1
  var ca1 = ee.Array(cm.consumersAccuracy()).get([0, 1]);  // precision class-1
  var f1  = pa1.multiply(ca1).multiply(2).divide(pa1.add(ca1));

  print('─────────────────────────────────────');
  print(label);
  print('Confusion matrix:',   cm);
  print('Overall accuracy:',   cm.accuracy());
  print('Kappa:',              cm.kappa());
  print('Producers accuracy:', cm.producersAccuracy());
  print('Consumers accuracy:', cm.consumersAccuracy());
  print('F1 (class retama):',  f1);
};

// =============================================================================
// 3. SECTION B — EXTERNAL VALIDATION 2023
// -----------------------------------------------------------------------------
// sampleRegions samples each prediction raster at the validation point locations.
// Points that fall on nodata pixels (e.g. RF2 after patch masking) are dropped
// automatically (dropNulls default = true). The count of matched points is
// printed so the drop rate can be assessed.
// =============================================================================
print('═════════════════════════════════════');
print('SECTION B — External validation ' + VAL_YEAR);

// --- B1. RF1 ---
var rf1_sampled = rf1_pred.sampleRegions({
  collection : valPoints,
  properties : ['label'],
  scale      : 10,
  tileScale  : 4
});
print('B1 RF1 — matched points:', rf1_sampled.size());

var cm_rf1_ext = rf1_sampled.errorMatrix({
  actual   : 'label',
  predicted: 'pred',
  order    : [0, 1]
});
printMetrics('RF1 — External validation ' + VAL_YEAR, cm_rf1_ext);

// --- B2. RF2 raw (before patch filter) ---
var rf2_raw_sampled = rf2_raw_pred.sampleRegions({
  collection : valPoints,
  properties : ['label'],
  scale      : 10,
  tileScale  : 4
});
print('B2 RF2raw — matched points:', rf2_raw_sampled.size());

var cm_rf2_raw_ext = rf2_raw_sampled.errorMatrix({
  actual   : 'label',
  predicted: 'classification',
  order    : [0, 1]
});
printMetrics('RF2 raw — External validation ' + VAL_YEAR, cm_rf2_raw_ext);

// --- B3. RF2 final (after patch filter) ---
var rf2_sampled = rf2_pred.sampleRegions({
  collection : valPoints,
  properties : ['label'],
  scale      : 10,
  tileScale  : 4
});
print('B3 RF2 final — matched points:', rf2_sampled.size());

var cm_rf2_ext = rf2_sampled.errorMatrix({
  actual   : 'label',
  predicted: 'classification',
  order    : [0, 1]
});
printMetrics('RF2 final — External validation ' + VAL_YEAR, cm_rf2_ext);

// =============================================================================
// 4. EXPORT SAMPLED TABLES TO DRIVE (for R analysis)
// -----------------------------------------------------------------------------
// Each exported CSV contains: label (0/1), predicted value, TIPO, YEAR.
// Allows recomputing any metric in R (e.g. ROC, precision-recall curves).
// =============================================================================
Export.table.toDrive({
  collection     : rf1_sampled,
  description    : 'Validation_RF1_external_' + VAL_YEAR,
  folder         : DRIVE_FOLDER,
  fileNamePrefix : 'Validation_RF1_external_' + VAL_YEAR,
  fileFormat     : 'CSV'
});

Export.table.toDrive({
  collection     : rf2_raw_sampled,
  description    : 'Validation_RF2raw_external_' + VAL_YEAR + RUN_SUFFIX,
  folder         : DRIVE_FOLDER,
  fileNamePrefix : 'Validation_RF2raw_external_' + VAL_YEAR + RUN_SUFFIX,
  fileFormat     : 'CSV'
});

Export.table.toDrive({
  collection     : rf2_sampled,
  description    : 'Validation_RF2final_external_' + VAL_YEAR + RUN_SUFFIX,
  folder         : DRIVE_FOLDER,
  fileNamePrefix : 'Validation_RF2final_external_' + VAL_YEAR + RUN_SUFFIX,
  fileFormat     : 'CSV'
});

// =============================================================================
// 5. SECTION C — CROSS-RUN COMPARISON
// -----------------------------------------------------------------------------
// Loops over COMPARE_SUFFIXES.  For each run:
//   C-B2 — RF2 raw   external validation VAL_YEAR
//   C-B3 — RF2 final external validation VAL_YEAR
// RF1 is run-independent → already shown in Section B, not repeated here.
// Console output only (no exports — use sections B-4 for the primary run).
// =============================================================================
print('═════════════════════════════════════');
print('SECTION C — Cross-run comparison (external ' + VAL_YEAR + ')');
print('Runs compared:', COMPARE_SUFFIXES.length);

COMPARE_SUFFIXES.forEach(function(suffix) {
  print('');
  print('▶▶▶ Run: ' + suffix);

  // --- C-B2. RF2 raw external ---
  var cmp_rf2_raw = ee.Image(ASSET_PREFIX + 'RF2_raw_prediction_' + VAL_YEAR + suffix)
    .select('classification');

  var cmp_raw_sampled = cmp_rf2_raw.sampleRegions({
    collection : valPoints,
    properties : ['label'],
    scale      : 10,
    tileScale  : 4
  });
  print('C-B2 [' + suffix + '] RF2 raw — matched points:', cmp_raw_sampled.size());

  var cmp_cm_raw = cmp_raw_sampled.errorMatrix({
    actual   : 'label',
    predicted: 'classification',
    order    : [0, 1]
  });
  printMetrics('[' + suffix + '] RF2 raw — External validation ' + VAL_YEAR, cmp_cm_raw);

  // --- C-B3. RF2 final external ---
  var cmp_rf2_fin = ee.Image(ASSET_PREFIX + 'RF2_prediction_' + VAL_YEAR + suffix)
    .select('classification');

  var cmp_fin_sampled = cmp_rf2_fin.sampleRegions({
    collection : valPoints,
    properties : ['label'],
    scale      : 10,
    tileScale  : 4
  });
  print('C-B3 [' + suffix + '] RF2 final — matched points:', cmp_fin_sampled.size());

  var cmp_cm_fin = cmp_fin_sampled.errorMatrix({
    actual   : 'label',
    predicted: 'classification',
    order    : [0, 1]
  });
  printMetrics('[' + suffix + '] RF2 final — External validation ' + VAL_YEAR, cmp_cm_fin);
});
