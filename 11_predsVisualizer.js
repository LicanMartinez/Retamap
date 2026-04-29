/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var trueColorVisParam = {"opacity":1,"bands":["B4","B3","B2"],"min":365.91480931730996,"max":1564.904287792141,"gamma":0.702};
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// =============================================================================
// 11: Predictions Visualizer — RF1 / RF2 raw / RF2 gap-fill / RF2 final, all years
// -----------------------------------------------------------------------------
// Prerequisite : assets 02_RF1_raw_prediction_YYYY        (from 02_RF1fit.js)
//              + assets 06_RF2_raw_prediction_YYYY{SUFFIX} (from 06_RF2predict.js)
//              + assets 07_RF2_gapFill_YYYY{SUFFIX}        (from 07_RF2gapFill.js)
//              + assets 08_RF2_prediction_YYYY{SUFFIX}     (from 08_RF2patchFilter.js)
// -----------------------------------------------------------------------------
// Displays only class-1 pixels (selfMask) for each model and year:
//   RF1          → yellow     (#FFD700)
//   RF2 raw      → orange     (#FF8C00)
//   RF2 gap-fill → teal       (#00BCD4)
//   RF2 final    → red        (#CC0000)
//
// Each layer is off by default; toggle in the Layers panel.
// =============================================================================

// =============================================================================
// 0. CONFIGURATION — edit years here
// =============================================================================
var ASSET_PREFIX  = 'projects/ee-licanemartinez/assets/Retamap/';
// var RUN_SUFFIX    = '_trimmedRF1comp_s1n12k_qs1n12k_s0n12k_qs0n12k';
var RUN_SUFFIX    = '_trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k';
var years      = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
// var years         = [2023, 2018];
var roi           = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');

// Set to true only after running the corresponding script
var SHOW_RF2_GAPFILL = true;   // toggle after 07_RF2gapFill.js
var SHOW_RF2_FINAL   = true;   // toggle after 08_RF2patchFilter.js

// Visualization params — only class-1 shown (selfMask removes class-0)
var visRF1     = {min: 1, max: 1, palette: ['#FFD700']};  // yellow
var visRF1_compCtrl     = {min: 1, max: 1, palette: ['green']};  // yellow
var visRF1_raw     = {min: 1, max: 1, palette: ['violet']};  // violet
var visRF2raw  = {min: 1, max: 1, palette: ['#FF8C00']};  // orange
var visRF2gf   = {min: 1, max: 1, palette: ['#00BCD4']};  // teal
var visRF2fin  = {min: 1, max: 1, palette: ['#CC0000']};  // red

// =============================================================================
// 2. ADD LAYERS
// -----------------------------------------------------------------------------
// Order per year (bottom → top): mosaic true color → NDYI → RF1 → RF2 raw
//                                 → RF2 gap-fill → RF2 final.
// All layers off by default.
// =============================================================================
// Map.centerObject(roi, 10);

// Mosaic visualization params (01_MergedBands_YYYY bands)
// var visTrueColor = {bands: ['B4', 'B3', 'B2'], min: 0, max: 2500};
var visTrueColor = trueColorVisParam
var visNDYI      = {min: -0.1, max: 0.4, palette: ['#313695','#abd9e9','#ffffbf','#fdae61','#f46d43','#d73027']};

years.forEach(function(year) {
  var mosaic = ee.Image(ASSET_PREFIX + '01_MergedBands_' + year).clip(roi);

  Map.addLayer(mosaic, visTrueColor, 'True color ' + year, false);
  Map.addLayer(mosaic.select('NDYI'), visNDYI, 'NDYI ' + year, false);

  // RF1 band is 'pred' (from 02_RF1fit.js)
  var rf1 = ee.Image(ASSET_PREFIX + '02_RF1_raw_prediction_' + year)
    .select('pred').selfMask().clip(roi);

  var rf1_compCtrl = ee.Image(ASSET_PREFIX + '02_RF1_2compCTRL_prediction_connectedFilter_' + year)
    .select('pred').selfMask().clip(roi);

  // RF2 raw band is 'classification' (from 06_RF2predict.js)
  var rf2raw = ee.Image(ASSET_PREFIX + '06_RF2_raw_prediction_' + year + RUN_SUFFIX)
    .select('classification').selfMask().clip(roi);

  // RF2 gap-fill band is 'classification' (from 07_RF2gapFill.js)
  var rf2gf = ee.Image(ASSET_PREFIX + '07_RF2_gapFill_' + year + RUN_SUFFIX)
    .select('classification').selfMask().clip(roi);

  // RF2 final band is 'classification' (from 08_RF2patchFilter.js)
  var rf2fin = ee.Image(ASSET_PREFIX + '08_RF2_prediction_' + year + RUN_SUFFIX)
    .select('classification').selfMask().clip(roi);

  Map.addLayer(rf1,    visRF1,    'RF1 '         + year, false);
  Map.addLayer(rf1_compCtrl,    visRF1_compCtrl,    'RF1_compCtrl '         + year, false);
  Map.addLayer(rf2raw, visRF2raw, 'RF2 raw '     + year, false);
  if (SHOW_RF2_GAPFILL) {
    Map.addLayer(rf2gf, visRF2gf, 'RF2 gap-fill ' + year, false);
  }
  if (SHOW_RF2_FINAL) {
    Map.addLayer(rf2fin, visRF2fin, 'RF2 final '  + year, false);
  }
});


// =============================================================================
// 3. MULTI-RUN RAW COMPARISON
// -----------------------------------------------------------------------------
// Compares RF2 raw predictions across different RUN_SUFFIX values + RF1.
// Edit COMPARE_SUFFIXES to include the runs you want to compare.
// Each run gets a distinct color; RF1 always shown in yellow.
// All layers off by default.
// =============================================================================

// var COMPARE_YEARS    = [2017, 2018];
var COMPARE_YEARS    = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var COMPARE_SUFFIXES = [
  '_trimmedRF1comp_s1n12k_qs1n12k_s0n12k_qs0n12k',
];

// Palettes for each suffix — extend if adding more runs
var COMPARE_PALETTES = [
  '#FF8C00',   // orange (run 1)
  '#00BFFF',   // deep-sky-blue (run 2)
  '#39FF14',   // neon-green (run 3)
  '#FF00FF',   // magenta (run 4)
  '#FF4500',   // orange-red (run 5)
];

// COMPARE_YEARS.forEach(function(year) {
//   // RF1 — shown once per year regardless of suffix list
//   // var rf1 = ee.Image(ASSET_PREFIX + '02_RF1_prediction_connectedFilter_' + year)
//   //   .select('pred').selfMask().clip(roi);
//   // Map.addLayer(rf1, visRF1, '[CMP] RF1 conn' + year, false);
  
//   // var rf1_raw = ee.Image(ASSET_PREFIX + '02_RF1_raw_prediction_' + year)
//   //   .select('pred').selfMask().clip(roi);
//   // Map.addLayer(rf1_raw, visRF1_raw, '[CMP] RF1 raw ' + year, false);
  
//   var rf1_comp = ee.Image(ASSET_PREFIX + '02_RF1_2compCTRL_prediction_connectedFilter_' + year)
//     .select('pred').selfMask().clip(roi);
//   Map.addLayer(rf1_comp, {min: 1, max: 1, palette: ['red']}, '[CMP] RF1 compCTRL ' + year, false);

//   // RF2 raw per suffix
//   COMPARE_SUFFIXES.forEach(function(suffix, idx) {
//     var palette = COMPARE_PALETTES[idx % COMPARE_PALETTES.length];
//     var rf2raw = ee.Image(ASSET_PREFIX + '06_RF2_raw_prediction_' + year + suffix)
//       .select('classification').selfMask().clip(roi);
//     Map.addLayer(
//       rf2raw,
//       {min: 1, max: 1, palette: [palette]},
//       '[CMP] RF2 raw ' + year + ' ' + suffix,
//       false
//     );
//   });
// });

// =============================================================================

// Print validation points
var valPoints = ee.FeatureCollection(ASSET_PREFIX + '4-Validation_points_complete_2023');

var ctrlPoints   = valPoints.filter(ee.Filter.eq('TIPO', 'ctrl'));
var retamaPoints = valPoints.filter(ee.Filter.eq('TIPO', 'retama'));

Map.addLayer(ctrlPoints,   {color: 'darkgreen'}, 'Puntos Ctrl', false);
Map.addLayer(retamaPoints, {color: 'gold'},       'Puntos Retama', false);
Map.setOptions('SATELLITE');
