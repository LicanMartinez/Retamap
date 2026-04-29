// =============================================================================
// 10: Validation — external validation 2023
// -----------------------------------------------------------------------------
// Prerequisite : asset 02_RF1_raw_prediction_2023          (from 02_RF1fit.js)
//              + asset 06_RF2_raw_prediction_2023{SUFFIX}  (from 06_RF2predict.js)
//              + asset 07_RF2_gapFill_2023{SUFFIX}         (from 07_RF2gapFill.js)
//              + asset 08_RF2_prediction_2023{SUFFIX}      (from 08_RF2patchFilter.js)
//              + asset 4-Validation_points_complete_2023
// Produces     : console output with full metrics for all stages
//              + Export.table.toDrive (4 CSV) for downstream analysis in R
// -----------------------------------------------------------------------------
// Validation structure
// └─ B. External validation 2023 (independent points)
//       ├─ B1. RF1 raw (02_RF1_raw_prediction_2023)
//       ├─ B2. RF2 raw (06_RF2_raw_prediction_2023)
//       ├─ B3. RF2 gap-fill (07_RF2_gapFill_2023)
//       └─ B4. RF2 final (08_RF2_prediction_2023)
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER    = 'Retamap/Retamap_GEE_Exports';
var EXPORT_TO_DRIVE = false;  // toggle: true → also export to Google Drive

var VAL_YEAR   = 2023;
var RUN_SUFFIX = '_trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k';

// Section C — cross-run comparison (edit to include the runs you want to compare)
var COMPARE_SUFFIXES = [
  '_s1.10k_qs1.0_s0.20k_qs0.0',
];

// =============================================================================
// 1. LOAD ASSETS
// =============================================================================

// RF1 uses band 'pred'; RF2 scripts use band 'classification'
var rf1_pred     = ee.Image(ASSET_PREFIX + '02_RF1_raw_prediction_' + VAL_YEAR).select('pred');
var rf2_raw_pred = ee.Image(ASSET_PREFIX + '06_RF2_raw_prediction_' + VAL_YEAR + RUN_SUFFIX).select('classification');
var rf2_gf_pred  = ee.Image(ASSET_PREFIX + '07_RF2_gapFill_'        + VAL_YEAR + RUN_SUFFIX).select('classification');
var rf2_pred     = ee.Image(ASSET_PREFIX + '08_RF2_prediction_'      + VAL_YEAR + RUN_SUFFIX).select('classification');

// External validation points
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
// 3. SECTION B — EXTERNAL VALIDATION 2023
// =============================================================================
print('═════════════════════════════════════');
print('SECTION B — External validation ' + VAL_YEAR);

// --- B1. RF1 ---
var rf1_sampled = rf1_pred.sampleRegions({
  collection: valPoints,
  properties: ['label'],
  scale     : 10,
  tileScale : 4
});
print('B1 RF1 — matched points:', rf1_sampled.size());

var cm_rf1_ext = rf1_sampled.errorMatrix({
  actual   : 'label',
  predicted: 'pred',
  order    : [0, 1]
});

// --- B2. RF2 raw ---
var rf2_raw_sampled = rf2_raw_pred.sampleRegions({
  collection: valPoints,
  properties: ['label'],
  scale     : 10,
  tileScale : 4
});
print('B2 RF2 raw — matched points:', rf2_raw_sampled.size());

var cm_rf2_raw_ext = rf2_raw_sampled.errorMatrix({
  actual   : 'label',
  predicted: 'classification',
  order    : [0, 1]
});

// --- B3. RF2 gap-fill ---
var rf2_gf_sampled = rf2_gf_pred.sampleRegions({
  collection: valPoints,
  properties: ['label'],
  scale     : 10,
  tileScale : 4
});
print('B3 RF2 gap-fill — matched points:', rf2_gf_sampled.size());

var cm_rf2_gf_ext = rf2_gf_sampled.errorMatrix({
  actual   : 'label',
  predicted: 'classification',
  order    : [0, 1]
});

// --- B4. RF2 final (after patch filter) ---
var rf2_sampled = rf2_pred.sampleRegions({
  collection: valPoints,
  properties: ['label'],
  scale     : 10,
  tileScale : 4
});
print('B4 RF2 final — matched points:', rf2_sampled.size());

var cm_rf2_ext = rf2_sampled.errorMatrix({
  actual   : 'label',
  predicted: 'classification',
  order    : [0, 1]
});

// --- B. Summary table ---
var sectionB_FC = ee.FeatureCollection([
  metricsFeature('RF1',          cm_rf1_ext),
  metricsFeature('RF2 raw',      cm_rf2_raw_ext),
  metricsFeature('RF2 gap-fill', cm_rf2_gf_ext),
  metricsFeature('RF2 final',    cm_rf2_ext)
]);

print(
  ui.Chart.feature.byFeature(sectionB_FC, 'Model', ['OA', 'Kappa', 'Recall', 'Precision', 'F1'])
    .setChartType('Table')
    .setOptions({title: 'Section B — External validation ' + VAL_YEAR})
);

// =============================================================================
// 4. EXPORT SAMPLED TABLES TO DRIVE (for R analysis)
// =============================================================================
if (EXPORT_TO_DRIVE) {
  Export.table.toDrive({
    collection    : rf1_sampled,
    description   : 'Validation_RF1_external_' + VAL_YEAR,
    folder        : DRIVE_FOLDER,
    fileNamePrefix: 'Validation_RF1_external_' + VAL_YEAR,
    fileFormat    : 'CSV'
  });

  Export.table.toDrive({
    collection    : rf2_raw_sampled,
    description   : 'Validation_RF2raw_external_' + VAL_YEAR + RUN_SUFFIX,
    folder        : DRIVE_FOLDER,
    fileNamePrefix: 'Validation_RF2raw_external_' + VAL_YEAR + RUN_SUFFIX,
    fileFormat    : 'CSV'
  });

  Export.table.toDrive({
    collection    : rf2_gf_sampled,
    description   : 'Validation_RF2gapFill_external_' + VAL_YEAR + RUN_SUFFIX,
    folder        : DRIVE_FOLDER,
    fileNamePrefix: 'Validation_RF2gapFill_external_' + VAL_YEAR + RUN_SUFFIX,
    fileFormat    : 'CSV'
  });

  Export.table.toDrive({
    collection    : rf2_sampled,
    description   : 'Validation_RF2final_external_' + VAL_YEAR + RUN_SUFFIX,
    folder        : DRIVE_FOLDER,
    fileNamePrefix: 'Validation_RF2final_external_' + VAL_YEAR + RUN_SUFFIX,
    fileFormat    : 'CSV'
  });
}

// =============================================================================
// 5. SECTION C — CROSS-RUN COMPARISON
// -----------------------------------------------------------------------------
// For each run: RF2 raw (C-B2), RF2 gap-fill (C-B3), RF2 final (C-B4).
// RF1 is run-independent → shown in Section B, not repeated here.
// =============================================================================
print('═════════════════════════════════════');
print('SECTION C — Cross-run comparison (external ' + VAL_YEAR + ')');
print('Runs compared:', COMPARE_SUFFIXES.length);

var sectionC_features = [];

COMPARE_SUFFIXES.forEach(function(suffix) {
  // --- C-B2. RF2 raw ---
  var cmp_raw = ee.Image(ASSET_PREFIX + '06_RF2_raw_prediction_' + VAL_YEAR + suffix)
    .select('classification');
  var cmp_raw_sampled = cmp_raw.sampleRegions({
    collection: valPoints, properties: ['label'], scale: 10, tileScale: 4
  });
  print('C-B2 [' + suffix + '] RF2 raw — matched:', cmp_raw_sampled.size());
  var cmp_cm_raw = cmp_raw_sampled.errorMatrix({
    actual: 'label', predicted: 'classification', order: [0, 1]
  });
  sectionC_features.push(metricsFeature(suffix + ' · raw', cmp_cm_raw));

  // --- C-B3. RF2 gap-fill ---
  var cmp_gf = ee.Image(ASSET_PREFIX + '07_RF2_gapFill_' + VAL_YEAR + suffix)
    .select('classification');
  var cmp_gf_sampled = cmp_gf.sampleRegions({
    collection: valPoints, properties: ['label'], scale: 10, tileScale: 4
  });
  print('C-B3 [' + suffix + '] RF2 gap-fill — matched:', cmp_gf_sampled.size());
  var cmp_cm_gf = cmp_gf_sampled.errorMatrix({
    actual: 'label', predicted: 'classification', order: [0, 1]
  });
  sectionC_features.push(metricsFeature(suffix + ' · gap-fill', cmp_cm_gf));

  // --- C-B4. RF2 final ---
  var cmp_fin = ee.Image(ASSET_PREFIX + '08_RF2_prediction_' + VAL_YEAR + suffix)
    .select('classification');
  var cmp_fin_sampled = cmp_fin.sampleRegions({
    collection: valPoints, properties: ['label'], scale: 10, tileScale: 4
  });
  print('C-B4 [' + suffix + '] RF2 final — matched:', cmp_fin_sampled.size());
  var cmp_cm_fin = cmp_fin_sampled.errorMatrix({
    actual: 'label', predicted: 'classification', order: [0, 1]
  });
  sectionC_features.push(metricsFeature(suffix + ' · final', cmp_cm_fin));
});

print(
  ui.Chart.feature.byFeature(
    ee.FeatureCollection(sectionC_features), 'Model', ['OA', 'Kappa', 'Recall', 'Precision', 'F1']
  )
    .setChartType('Table')
    .setOptions({title: 'Section C — Cross-run comparison (external ' + VAL_YEAR + ')'})
);
