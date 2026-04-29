// =============================================================================
// 09: Holdout validation — evaluate RF2 stages against held-out val samples
// -----------------------------------------------------------------------------
// Prerequisite : assets 04_RF2_samples_YYYY{SUFFIX}_val  (from 04_RF2sampleExport.js)
//              + assets 06_RF2_raw_prediction_YYYY{SUFFIX}  (from 06_RF2predict.js)
//              + assets 07_RF2_gapFill_YYYY{SUFFIX}         (from 07_RF2gapFill.js)
//              + assets 08_RF2_prediction_YYYY{SUFFIX}      (from 08_RF2patchFilter.js)
// Produces     : console summary tables per stage
//              + optional CSV exports to Drive (EXPORT_TO_DRIVE = true)
// -----------------------------------------------------------------------------
// Validation structure (MOD 5):
//   For each year in VAL_YEARS:
//     Stage A  — RF2 raw          (06_RF2_raw_prediction_YYYY)
//     Stage B  — RF2 gap-fill     (07_RF2_gapFill_YYYY)
//     Stage C  — RF2 final        (08_RF2_prediction_YYYY)
//   Summary table aggregated across all years is printed at the end.
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER    = 'Retamap/Retamap_GEE_Exports';
var EXPORT_TO_DRIVE = false;  // toggle: true → export matched-point tables to Drive

var exportYears = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var VAL_YEARS   = exportYears;   // edit to evaluate a subset only
var RUN_SUFFIX  = '_trimmedRF1comp_s1n12k_qs1n12k_s0n12k_qs0n12k';

// =============================================================================
// 1. HELPER: metrics feature from confusion matrix
// =============================================================================
var metricsFeature = function(label, cm, count) {
  var pa1 = ee.Array(cm.producersAccuracy()).get([1, 0]);  // recall    class-1
  var ca1 = ee.Array(cm.consumersAccuracy()).get([0, 1]);  // precision class-1
  var f1  = pa1.multiply(ca1).multiply(2).divide(pa1.add(ca1));
  return ee.Feature(null, {
    'Stage'         : label,
    'Matched_Points': count,
    'OA'            : cm.accuracy(),
    'Kappa'         : cm.kappa(),
    'Recall'        : pa1,
    'Precision'     : ca1,
    'F1'            : f1
  });
};

// =============================================================================
// 2. PER-YEAR VALIDATION
// =============================================================================
print('══════════════════════════════════════════');
print('09 — Holdout validation (' + RUN_SUFFIX + ')');

var allFeatures = [];

VAL_YEARS.forEach(function(year) {
  print('────── Year ' + year + ' ──────');

  // Load val samples for this year
  var valSamples = ee.FeatureCollection(
    ASSET_PREFIX + '04_RF2_samples_' + year + RUN_SUFFIX + '_val'
  );
  print('  Val samples (' + year + '):', valSamples.size());

  // ── Stage A: RF2 raw ────────────────────────────────────────────────────────
  var raw = ee.Image(ASSET_PREFIX + '06_RF2_raw_prediction_' + year + RUN_SUFFIX)
    .select('classification');

  var raw_sampled = raw.sampleRegions({
    collection: valSamples,
    properties: ['stable_label'],
    scale     : 10,
    tileScale : 4
  });
  print('  A RF2 raw — matched:', raw_sampled.size());

  var cm_raw = raw_sampled.errorMatrix({
    actual   : 'stable_label',
    predicted: 'classification',
    order    : [0, 1]
  });
  allFeatures.push(metricsFeature(year + ' A raw', cm_raw, raw_sampled.size()));

  if (EXPORT_TO_DRIVE) {
    Export.table.toDrive({
      collection    : raw_sampled,
      description   : 'Holdout_RF2raw_' + year + RUN_SUFFIX,
      folder        : DRIVE_FOLDER,
      fileNamePrefix: 'Holdout_RF2raw_' + year + RUN_SUFFIX,
      fileFormat    : 'CSV'
    });
  }

  // ── Stage B: RF2 gap-fill ───────────────────────────────────────────────────
  var gf = ee.Image(ASSET_PREFIX + '07_RF2_gapFill_' + year + RUN_SUFFIX)
    .select('classification');

  var gf_sampled = gf.sampleRegions({
    collection: valSamples,
    properties: ['stable_label'],
    scale     : 10,
    tileScale : 4
  });
  print('  B RF2 gap-fill — matched:', gf_sampled.size());

  var cm_gf = gf_sampled.errorMatrix({
    actual   : 'stable_label',
    predicted: 'classification',
    order    : [0, 1]
  });
  allFeatures.push(metricsFeature(year + ' B gap-fill', cm_gf, gf_sampled.size()));

  if (EXPORT_TO_DRIVE) {
    Export.table.toDrive({
      collection    : gf_sampled,
      description   : 'Holdout_RF2gapFill_' + year + RUN_SUFFIX,
      folder        : DRIVE_FOLDER,
      fileNamePrefix: 'Holdout_RF2gapFill_' + year + RUN_SUFFIX,
      fileFormat    : 'CSV'
    });
  }

  // ── Stage C: RF2 final (patch-filtered) ────────────────────────────────────
  var fin = ee.Image(ASSET_PREFIX + '08_RF2_prediction_' + year + RUN_SUFFIX)
    .select('classification');

  var fin_sampled = fin.sampleRegions({
    collection: valSamples,
    properties: ['stable_label'],
    scale     : 10,
    tileScale : 4
  });
  print('  C RF2 final — matched:', fin_sampled.size());

  var cm_fin = fin_sampled.errorMatrix({
    actual   : 'stable_label',
    predicted: 'classification',
    order    : [0, 1]
  });
  allFeatures.push(metricsFeature(year + ' C final', cm_fin, fin_sampled.size()));

  if (EXPORT_TO_DRIVE) {
    Export.table.toDrive({
      collection    : fin_sampled,
      description   : 'Holdout_RF2final_' + year + RUN_SUFFIX,
      folder        : DRIVE_FOLDER,
      fileNamePrefix: 'Holdout_RF2final_' + year + RUN_SUFFIX,
      fileFormat    : 'CSV'
    });
  }
});

// =============================================================================
// 3. SUMMARY TABLE (all years × all stages)
// =============================================================================
print('══════════════════════════════════════════');
print('Summary — all years and stages');
print(
  ui.Chart.feature.byFeature(
    ee.FeatureCollection(allFeatures),
    'Stage',
    ['Matched_Points', 'OA', 'Kappa', 'Recall', 'Precision', 'F1']
  )
    .setChartType('Table')
    .setOptions({title: '09 — Holdout validation (' + RUN_SUFFIX + ')'})
);
