// =============================================================================
// 09.2: Holdout validation — standard remote-sensing metrics
// -----------------------------------------------------------------------------
// Prerequisite : assets 04_RF2_samples_YYYY{SUFFIX}_val  (from 04_RF2sampleExport.js)
//              + assets 06_RF2_raw_prediction_YYYY{SUFFIX}  (from 06_RF2predict.js)
//              + assets 07_RF2_gapFill_YYYY{SUFFIX}         (from 07_RF2gapFill.js)
//              + assets 08_RF2_prediction_YYYY{SUFFIX}      (from 08_RF2patchFilter.js)
// Produces     : console summary table (all years × all stages)
//              + optional CSV exports to Drive (EXPORT_TO_DRIVE = true)
// -----------------------------------------------------------------------------
// Metrics reported (both classes, as per RS standard):
//   OA, Kappa, PA_retama, UA_retama, OE_retama, CE_retama,
//   PA_ctrl, UA_ctrl, F1_retama, TP, FP, FN, TN
//
// Stage C (patch-filtered) note:
//   Pixels removed by the patch filter are class-0 predictions (the filter
//   decided they are not retama). They are recovered via unmask(0) within the
//   gap-fill valid mask so all stages are evaluated on the same N points.
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER    = 'Retamap/Retamap_GEE_Exports';
var EXPORT_TO_DRIVE = false;

var exportYears = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var VAL_YEARS   = exportYears;
var RUN_SUFFIX  = '_trimmedRF1comp_s1n12k_qs1n12k_s0n12k_qs0n12k';
// var RUN_SUFFIX  = '_trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k';

// =============================================================================
// 1. HELPER: full RS metrics from confusion matrix
// =============================================================================
// Confusion matrix layout with order=[0,1]:  [[TN, FP], [FN, TP]]
var metricsFeature = function(label, cm, n) {
  var mat = cm.array();
  var TN  = mat.get([0, 0]);
  var FP  = mat.get([0, 1]);
  var FN  = mat.get([1, 0]);
  var TP  = mat.get([1, 1]);

  var pa1 = TP.divide(TP.add(FN));          // Producer's Accuracy class-1 (Recall)
  var ua1 = TP.divide(TP.add(FP));          // User's Accuracy    class-1 (Precision)
  var pa0 = TN.divide(TN.add(FP));          // Producer's Accuracy class-0
  var ua0 = TN.divide(TN.add(FN));          // User's Accuracy    class-0
  var f1  = pa1.multiply(ua1).multiply(2).divide(pa1.add(ua1));

  return ee.Feature(null, {
    'Stage'    : label,
    'N'        : n,
    'OA'       : cm.accuracy(),
    'Kappa'    : cm.kappa(),
    'PA_retama': pa1,
    'UA_retama': ua1,
    'OE_retama': ee.Number(1).subtract(pa1),  // Omission Error  class-1
    'CE_retama': ee.Number(1).subtract(ua1),  // Commission Error class-1
    'PA_ctrl'  : pa0,
    'UA_ctrl'  : ua0,
    'F1'       : f1,
    'TP': TP, 'FP': FP, 'FN': FN, 'TN': TN
  });
};

var sampleAndMetrics = function(image, valSamples, predBand, label) {
  var sampled = image.sampleRegions({
    collection: valSamples,
    properties: ['stable_label'],
    scale     : 10,
    tileScale : 4
  });
  var cm = sampled.errorMatrix({
    actual   : 'stable_label',
    predicted: predBand,
    order    : [0, 1]
  });
  return metricsFeature(label, cm, sampled.size());
};

// =============================================================================
// 2. PER-YEAR VALIDATION
// =============================================================================
print('══════════════════════════════════════════════════');
print('09.2 — Holdout validation (' + RUN_SUFFIX + ')');

var allFeatures = [];

VAL_YEARS.forEach(function(year) {
  print('── ' + year + ' ──');

  var valSamples = ee.FeatureCollection(
    ASSET_PREFIX + '04_RF2_samples_' + year + RUN_SUFFIX + '_val'
  );

  var raw = ee.Image(ASSET_PREFIX + '06_RF2_raw_prediction_' + year + RUN_SUFFIX)
    .select('classification');

  var gf = ee.Image(ASSET_PREFIX + '07_RF2_gapFill_' + year + RUN_SUFFIX)
    .select('classification');

  // Stage C: recover patch-filtered pixels as class-0 within the valid mask.
  var fin = ee.Image(ASSET_PREFIX + '08_RF2_prediction_' + year + RUN_SUFFIX)
    .select('classification')
    .unmask(0)
    .updateMask(gf.mask());

  allFeatures.push(sampleAndMetrics(raw, valSamples, 'classification', year + ' A raw'));
  allFeatures.push(sampleAndMetrics(gf,  valSamples, 'classification', year + ' B gap-fill'));
  allFeatures.push(sampleAndMetrics(fin, valSamples, 'classification', year + ' C final'));

  if (EXPORT_TO_DRIVE) {
    var stages = [
      {img: raw, tag: 'RF2raw'},
      {img: gf,  tag: 'RF2gapFill'},
      {img: fin, tag: 'RF2final'}
    ];
    stages.forEach(function(s) {
      var sampled = s.img.sampleRegions({
        collection: valSamples, properties: ['stable_label'], scale: 10, tileScale: 4
      });
      Export.table.toDrive({
        collection    : sampled,
        description   : 'Holdout_' + s.tag + '_' + year + RUN_SUFFIX,
        folder        : DRIVE_FOLDER,
        fileNamePrefix: 'Holdout_' + s.tag + '_' + year + RUN_SUFFIX,
        fileFormat    : 'CSV'
      });
    });
  }
});

// =============================================================================
// 3. SUMMARY TABLE — all years × all stages
// =============================================================================
var cols = ['N', 'OA', 'Kappa',
            'PA_retama', 'UA_retama', 'OE_retama', 'CE_retama',
            'PA_ctrl', 'UA_ctrl', 'F1',
            'TP', 'FP', 'FN', 'TN'];

print('══════════════════════════════════════════════════');
print('Summary — all years and stages');
print(
  ui.Chart.feature.byFeature(ee.FeatureCollection(allFeatures), 'Stage', cols)
    .setChartType('Table')
    .setOptions({title: '09.2 — Holdout validation (' + RUN_SUFFIX + ')'})
);
