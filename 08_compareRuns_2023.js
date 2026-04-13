// =============================================================================
// 08_compareRuns_2023: Comparación de tres assets sobre validación externa 2023
// -----------------------------------------------------------------------------
// Compara directamente (sin holdout interno):
//   A. RF1_prediction_2023                        (band: pred)
//   B. RF2_raw_prediction_2023_noQS_n10k          (band: classification)
//   C. RF2_raw_prediction_2023_noQS1_n20k         (band: classification)
//
// Asset de validación: 4-Validation_points_complete_2023
//   Schema: TIPO ('retama'/'ctrl'), YEAR (2023)
// =============================================================================

// =============================================================================
// 0. ASSETS A COMPARAR
// =============================================================================
var ASSET_PREFIX = 'projects/ee-licanemartinez/assets/Retamap/';

var assetA = ee.Image(ASSET_PREFIX + 'RF1_prediction_2023')
               .select('pred').rename('classification');

var assetB = ee.Image(ASSET_PREFIX + 'RF2_raw_prediction_2023_noQS_n10k')
               .select('classification');

var assetC = ee.Image(ASSET_PREFIX + 'RF2_raw_prediction_2023_noQS1_n20k')
               .select('classification');

// =============================================================================
// 1. PUNTOS DE VALIDACIÓN
// =============================================================================
var valPointsRaw = ee.FeatureCollection(ASSET_PREFIX + '4-Validation_points_complete_2023');

// TIPO → label numérico (retama=1, ctrl=0)
var valPoints = valPointsRaw.map(function(f) {
  var label = ee.Algorithms.If(ee.String(f.get('TIPO')).equals('retama'), 1, 0);
  return f.set('label', label);
});

print('════════════════════════════════════');
print('Validation points (total):', valPoints.size());
print('  retama:', valPoints.filter(ee.Filter.eq('label', 1)).size());
print('  ctrl:  ', valPoints.filter(ee.Filter.eq('label', 0)).size());

// =============================================================================
// 2. HELPER: imprime métricas completas
// =============================================================================
var printMetrics = function(label, cm) {
  var pa1 = ee.Array(cm.producersAccuracy()).get([1, 0]);  // recall    clase-1
  var ca1 = ee.Array(cm.consumersAccuracy()).get([0, 1]);  // precision clase-1
  var f1  = pa1.multiply(ca1).multiply(2).divide(pa1.add(ca1));

  print('────────────────────────────────────');
  print(label);
  print('Matched points:',       ee.Number(cm.array().reduce(ee.Reducer.sum(), [0,1]).get([0,0])));
  print('Confusion matrix:',     cm);
  print('Overall accuracy:',     cm.accuracy());
  print('Kappa:',                cm.kappa());
  print('Producers accuracy:',   cm.producersAccuracy());
  print('Consumers accuracy:',   cm.consumersAccuracy());
  print('F1 (class retama):',    f1);
};

// =============================================================================
// 3. MUESTREO Y VALIDACIÓN
// =============================================================================
print('════════════════════════════════════');
print('External validation 2023');

// --- A. RF1 ---
var sampledA = assetA.sampleRegions({
  collection : valPoints,
  properties : ['label'],
  scale      : 10,
  tileScale  : 4
});
print('A — RF1_prediction_2023 — matched:', sampledA.size());
var cmA = sampledA.errorMatrix({actual: 'label', predicted: 'classification', order: [0, 1]});
printMetrics('A — RF1_prediction_2023', cmA);

// --- B. RF2 raw noQS n10k ---
var sampledB = assetB.sampleRegions({
  collection : valPoints,
  properties : ['label'],
  scale      : 10,
  tileScale  : 4
});
print('B — RF2_raw_prediction_2023_noQS_n10k — matched:', sampledB.size());
var cmB = sampledB.errorMatrix({actual: 'label', predicted: 'classification', order: [0, 1]});
printMetrics('B — RF2_raw_prediction_2023_noQS_n10k', cmB);

// --- C. RF2 raw noQS1 n20k ---
var sampledC = assetC.sampleRegions({
  collection : valPoints,
  properties : ['label'],
  scale      : 10,
  tileScale  : 4
});
print('C — RF2_raw_prediction_2023_noQS1_n20k — matched:', sampledC.size());
var cmC = sampledC.errorMatrix({actual: 'label', predicted: 'classification', order: [0, 1]});
printMetrics('C — RF2_raw_prediction_2023_noQS1_n20k', cmC);

// =============================================================================
// 4. EXPORT A DRIVE (opcional — descomentar para correr como tarea)
// =============================================================================
// Export.table.toDrive({collection: sampledA, description: 'Val_RF1_2023',              fileFormat: 'CSV'});
// Export.table.toDrive({collection: sampledB, description: 'Val_RF2raw_noQS_n10k_2023', fileFormat: 'CSV'});
// Export.table.toDrive({collection: sampledC, description: 'Val_RF2raw_noQS1_n20k_2023',fileFormat: 'CSV'});
