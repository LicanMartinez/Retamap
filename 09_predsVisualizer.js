// =============================================================================
// 09: Predictions Visualizer — RF1 / RF2 raw / RF2 final, all years
// -----------------------------------------------------------------------------
// Prerequisite : assets RF1_prediction_YYYY        (from 02_RF1fit.js)
//              + assets RF2_raw_prediction_YYYY     (from 06_RF2predict.js)
//              + assets RF2_prediction_YYYY         (from 07_RF2patchFilter.js)
// -----------------------------------------------------------------------------
// Displays only class-1 pixels (selfMask) for each model and year:
//   RF1        → yellow   (#FFD700)
//   RF2 raw    → orange   (#FF8C00)
//   RF2 final  → red      (#CC0000)
//
// Each layer is off by default; toggle in the Layers panel.
// =============================================================================

// =============================================================================
// 0. CONFIGURATION — edit years here
// =============================================================================
var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var RUN_SUFFIX      = '_QS_n10k';  // must match 06/07
// var years           = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var years           = [2023, 2017];
var roi             = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');

// Set to true only after running 07_RF2patchFilter.js
var SHOW_RF2_FINAL  = true;

// Visualization params — only class-1 shown (selfMask removes class-0)
var visRF1     = {min: 1, max: 1, palette: ['#FFD700']};  // yellow
var visRF2raw  = {min: 1, max: 1, palette: ['#FF8C00']};  // orange
var visRF2fin  = {min: 1, max: 1, palette: ['#CC0000']};  // red

// =============================================================================
// 2. ADD LAYERS
// -----------------------------------------------------------------------------
// Order per year: RF1 (bottom) → RF2 raw → RF2 final (top).
// All layers off by default.
// =============================================================================
// Map.centerObject(roi, 10);

years.forEach(function(year) {
  // RF1 band is 'pred' (from 02_RF1fit.js)
  var rf1 = ee.Image(ASSET_PREFIX + 'RF1_prediction_' + year)
    .select('pred').selfMask().clip(roi);

  // RF2 raw band is 'classification' (from 06_RF2predict.js)
  var rf2raw = ee.Image(ASSET_PREFIX + 'RF2_raw_prediction_' + year + RUN_SUFFIX)
    .select('classification').selfMask().clip(roi);

  // RF2 final band is 'classification' (from 07_RF2patchFilter.js)
  var rf2fin = ee.Image(ASSET_PREFIX + 'RF2_prediction_' + year + RUN_SUFFIX)
    .select('classification').selfMask().clip(roi);

  // Map.addLayer(rf1,    visRF1,    'RF1 '       + year, false);
  // Map.addLayer(rf2raw, visRF2raw, 'RF2 raw '   + year, false);
  // if (SHOW_RF2_FINAL) {
  //   Map.addLayer(rf2fin, visRF2fin, 'RF2 final ' + year, false);
  // }
});


// =============================================================================
// 3. MULTI-RUN RAW COMPARISON
// -----------------------------------------------------------------------------
// Compares RF2 raw predictions across different RUN_SUFFIX values + RF1.
// Edit COMPARE_SUFFIXES to include the runs you want to compare.
// Each run gets a distinct color; RF1 always shown in yellow.
// All layers off by default.
// =============================================================================

var COMPARE_YEARS    = [2023, 2017];            // years to compare (usually 1–2 key years)
var COMPARE_SUFFIXES = [                  // list of RUN_SUFFIX strings to compare
  '_noQS_n20k_compSamp',
  '_QS_n10k_compSamp',
  '_noQS1_n20k',
  '_noQS_n10k',
  '_QS_n10k'
  ];

// Palettes for each suffix — extend if adding more runs
var COMPARE_PALETTES = [
  '#FF8C00',   // orange (run 1)
  '#00BFFF',   // deep-sky-blue (run 2)
  '#39FF14',   // neon-green (run 3)
  '#FF00FF',   // magenta (run 4)
  '#FF4500',   // orange-red (run 5)
];

COMPARE_YEARS.forEach(function(year) {
  // RF1 — shown once per year regardless of suffix list
  var rf1 = ee.Image(ASSET_PREFIX + 'RF1_prediction_' + year)
    .select('pred').selfMask().clip(roi);
  Map.addLayer(rf1, visRF1, '[CMP] RF1 ' + year, false);

  // RF2 raw per suffix
  COMPARE_SUFFIXES.forEach(function(suffix, idx) {
    var palette = COMPARE_PALETTES[idx % COMPARE_PALETTES.length];
    var rf2raw = ee.Image(ASSET_PREFIX + 'RF2_raw_prediction_' + year + suffix)
      .select('classification').selfMask().clip(roi);
    Map.addLayer(
      rf2raw,
      {min: 1, max: 1, palette: [palette]},
      '[CMP] RF2 raw ' + year + ' ' + suffix,
      false
    );
  });
});

// =============================================================================

////// print validation points
// Cargar la FeatureCollection
var valPoints = ee.FeatureCollection('projects/ee-licanemartinez/assets/Retamap/5-Validation_points_complete_2023');

// Filtrar los puntos por el atributo 'TIPO'
var ctrlPoints = valPoints.filter(ee.Filter.eq('TIPO', 'ctrl'));
var retamaPoints = valPoints.filter(ee.Filter.eq('TIPO', 'retama'));

// Centrar el mapa en la extensión de los puntos
// Map.centerObject(valPoints, 10);

// Agregar las capas al mapa con los colores requeridos
Map.addLayer(ctrlPoints, {color: 'darkgreen'}, 'Puntos Ctrl');
Map.addLayer(retamaPoints, {color: 'gold'}, 'Puntos Retama');
Map.setOptions('SATELLITE')

