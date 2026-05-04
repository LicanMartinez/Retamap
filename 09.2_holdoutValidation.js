var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER    = 'Retamap';
var EXPORT_TO_DRIVE = true;

var exportYears = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var VAL_YEARS   = exportYears;
var RUN_SUFFIX  = '_trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k';

var metricsFeature = function(label, cm, n) {
  var mat = cm.array();
  var TN  = mat.get([0, 0]);
  var FP  = mat.get([0, 1]);
  var FN  = mat.get([1, 0]);
  var TP  = mat.get([1, 1]);

  var pa1 = TP.divide(TP.add(FN));
  var ua1 = TP.divide(TP.add(FP));
  var pa0 = TN.divide(TN.add(FP));
  var ua0 = TN.divide(TN.add(FN));
  var f1  = pa1.multiply(ua1).multiply(2).divide(pa1.add(ua1));

  return ee.Feature(null, {
    'Stage'    : label,
    'N'        : n,
    'OA'       : cm.accuracy(),
    'Kappa'    : cm.kappa(),
    'PA_retama': pa1,
    'UA_retama': ua1,
    'OE_retama': ee.Number(1).subtract(pa1),
    'CE_retama': ee.Number(1).subtract(ua1),
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

var allFeatures = [];

VAL_YEARS.forEach(function(year) {
  var valSamples = ee.FeatureCollection(
    ASSET_PREFIX + '04_RF2_samples_' + year + RUN_SUFFIX + '_val'
  );

  var raw = ee.Image(ASSET_PREFIX + '06_RF2_raw_prediction_' + year + RUN_SUFFIX)
    .select('classification');

  var gf = ee.Image(ASSET_PREFIX + '07_RF2_gapFill_' + year + RUN_SUFFIX)
    .select('classification');

  var fin = ee.Image(ASSET_PREFIX + '08_RF2_prediction_' + year + RUN_SUFFIX)
    .select('classification')
    .unmask(0)
    .updateMask(gf.mask());

  allFeatures.push(sampleAndMetrics(raw, valSamples, 'classification', year + ' A raw'));
  allFeatures.push(sampleAndMetrics(gf,  valSamples, 'classification', year + ' B gap-fill'));
  allFeatures.push(sampleAndMetrics(fin, valSamples, 'classification', year + ' C final'));
});

var metricsCollection = ee.FeatureCollection(allFeatures);

if (EXPORT_TO_DRIVE) {
  Export.table.toDrive({
    collection: metricsCollection,
    description: 'Holdout_Summary_Metrics' + RUN_SUFFIX,
    folder: DRIVE_FOLDER,
    fileNamePrefix: 'Holdout_Summary_Metrics' + RUN_SUFFIX,
    fileFormat: 'CSV'
  });
}