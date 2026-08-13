// =============================================================================
// 11.2: Predictions Visualizer — final map, full time series in one colour ramp
// -----------------------------------------------------------------------------
// Role         : TOOL — visualization only. Creates no asset, exports nothing.
// Prerequisite : assets 08_RF2_prediction_YYYY{SUFFIX} (from 08_RF2patchFilter.js)
// -----------------------------------------------------------------------------
// Years are drawn newest-first so the older ones sit on top: useful to see where
// the invasion gained and lost ground across the series.
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX  = 'projects/ee-licanemartinez/assets/Retamap/';
var RUN_SUFFIX    = '_trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k';
var roi           = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');

// Orden descendente: se agregan primero los nuevos. Los más antiguos se agregan 
// al final, por lo que quedan renderizados por encima de los demás.
var years = [2025, 2024, 2023, 2022, 2021, 2020, 2019, 2018, 2017];

// Colores: más oscuros para años antiguos, más claros para años recientes
var yearColors = {
  2017: '#800000', // Rojo oscuro (más antiguo)
  2018: '#b30000',
  2019: '#e60000', // Rojo
  2020: '#ff4400', // Rojo anaranjado
  2021: '#ff8800', // Naranja
  2022: '#ffbb00', // Naranja amarillento
  2023: '#ffdd00', // Dorado
  2024: '#ffff00', // Amarillo
  2025: '#ffff99'  // Amarillo claro (más reciente)
};

// =============================================================================
// 2. ADD LAYERS
// =============================================================================
// Map.centerObject(roi, 10);

years.forEach(function(year) {
  // RF2 final band is 'classification' (from 08_RF2patchFilter.js)
  var rf2fin = ee.Image(ASSET_PREFIX + '08_RF2_prediction_' + year + RUN_SUFFIX)
    .select('classification').selfMask().clip(roi);

  var visRF2fin = {min: 1, max: 1, palette: [yearColors[year]]};
  Map.addLayer(rf2fin, visRF2fin, 'RF2 final ' + year, true, 0.5);
});


// =============================================================================

Map.setOptions('SATELLITE');