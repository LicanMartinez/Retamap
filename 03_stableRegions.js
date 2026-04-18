// =============================================================================
// 03: Stable & quasi-stable regions from RF1 prediction assets
// Prerequisite : assets 02_RF1_raw_prediction_YYYY  (from 02_RF1fit.js)
// Produces     : asset 03_RF2_stable_categories
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER    = 'Retamap/Retamap_GEE_Exports';
var EXPORT_TO_DRIVE = false;  // toggle: true → also export to Google Drive

var exportYears = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var roi         = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');

// =============================================================================
// 1. LOAD RF1 PREDICTIONS AND STACK AS MULTI-BAND IMAGE
// =============================================================================
var predStack = ee.Image(
  exportYears.map(function(year) {
    return ee.Image(ASSET_PREFIX + '02_RF1_prediction_connectedFilter_' + year)
    // return ee.Image(ASSET_PREFIX + 'RF1_prediction_' + year)
      .select('pred')
      .rename('pred_' + year);
  })
);

var N = exportYears.length; // total number of years

// =============================================================================
// 2. SUM OF 1s ACROSS YEARS (per pixel)
// =============================================================================
var sumOnes = predStack.reduce(ee.Reducer.sum()).rename('sum_ones');


// =============================================================================
// 3. BUILD STABLE AND QUASI-STABLE MASKS
// =============================================================================

// Stable 1: all years == 1
var stable1      = sumOnes.eq(N);
// Quasi-stable 1: all years minus one == 1
var quasiStable1 = sumOnes.eq(N - 1);

// Stable 0: all years == 0
var stable0      = sumOnes.eq(0);
// Quasi-stable 0: all years minus one == 0 (i.e., sum == 1)
var quasiStable0 = sumOnes.eq(1);

// Combined mask: pixels that are stable or quasi-stable (either class)
var stableMask      = stable1.or(stable0);
var quasiStableMask = quasiStable1.or(quasiStable0);
var combinedMask    = stableMask.or(quasiStableMask);

// Labeled layer: 1=stable-retama, 2=quasi-stable-retama,
//                3=stable-non-retama, 4=quasi-stable-non-retama
var stableLayer = ee.Image(0)
  .where(stable0,      3)
  .where(quasiStable0, 4)
  .where(stable1,      1)
  .where(quasiStable1, 2)
  .updateMask(combinedMask)
  .rename('stable_class')
  .clip(roi);

// =============================================================================
// 4. HISTOGRAM OF sumOnes
// =============================================================================
var histDict = sumOnes.reduceRegion({
  reducer  : ee.Reducer.frequencyHistogram(),
  geometry : roi.geometry(),
  scale    : 10,
  maxPixels: 1e13,
  tileScale: 16
});

// Export CSV (03_SumOnes_Histogram)
var histValues = ee.Dictionary(histDict.get('sum_ones'));
var histFC = ee.FeatureCollection(
  histValues.keys().map(function(k) {
    return ee.Feature(null, {
      sum_ones   : ee.Number.parse(k),
      pixel_count: histValues.get(k)
    });
  })
);

if (EXPORT_TO_DRIVE) {
  Export.table.toDrive({
    collection    : histFC,
    description   : '03_SumOnes_Histogram_rf1patchFiltered',
    folder        : DRIVE_FOLDER,
    fileNamePrefix: '03_SumOnes_Histogram_rf1patchFiltered',
    fileFormat    : 'CSV'
  });
}

// Console chart
histDict.evaluate(function(d) {
  var hist = d['sum_ones'];
  var rows = Object.keys(hist)
    .map(Number)
    .sort(function(a, b) { return a - b; })
    .map(function(k) { return {sum_ones: k, pixels: hist[String(k)]}; });
  print('Histogram table:', rows);

  var chartData = rows.map(function(r) { return [r.sum_ones, r.pixels]; });
  chartData.unshift(['Sum of RF1=1', 'Pixel count']);

  var chart = ui.Chart(chartData)
    .setChartType('ColumnChart')
    .setOptions({
      title: 'Distribution of sumOnes across ROI',
      hAxis: {title: 'Years with RF1=1 (out of ' + N + ')', ticks: [0,1,2,3,4,5,6,7,8,9]},
      vAxis: {title: 'Pixel count', logScale: true},
      legend: {position: 'none'},
      colors: ['#2d9c00']
    });
  print(chart);
});

// =============================================================================
// 5. VISUALIZATION
// =============================================================================
// Map.centerObject(roi, 10);

Map.setOptions('SATELLITE')

Map.addLayer(sumOnes.selfMask().clip(roi),
  {min: 0, max: N, palette: ['white', 'yellow', 'orange', 'red', 'darkred']},
  'Sum of RF1=1 across years', false);

// Individual stable class layers
Map.addLayer(
  stable1.selfMask().clip(roi),
  {palette: ['#FFB300']},
  '1 - Stable retama'
);

Map.addLayer(
  quasiStable1.selfMask().clip(roi),
  {palette: ['#FFFF99']},
  '2 - Quasi-stable retama'
);

Map.addLayer(
  stable0.selfMask().clip(roi),
  {palette: ['#006400']},
  '3 - Stable non-retama'
);

Map.addLayer(
  quasiStable0.selfMask().clip(roi),
  {palette: ['#90EE90']},
  '4 - Quasi-stable non-retama'
);

// =============================================================================
// 6. INSPECTOR: anomaly year + click series
// =============================================================================

// --- 6a. Año atípico en píxeles quasi-estables ---

// Para qs1 (sumOnes == N-1): ¿en qué año único fue pred==0?
var qs1AnomalyYear = ee.Image(0);
exportYears.forEach(function(year) {
  qs1AnomalyYear = qs1AnomalyYear.where(
    predStack.select('pred_' + year).eq(0),
    year
  );
});
qs1AnomalyYear = qs1AnomalyYear
  .updateMask(quasiStable1)
  .rename('qs1_anomaly_year')
  .clip(roi);

// Para qs0 (sumOnes == 1): ¿en qué año único fue pred==1?
var qs0AnomalyYear = ee.Image(0);
exportYears.forEach(function(year) {
  qs0AnomalyYear = qs0AnomalyYear.where(
    predStack.select('pred_' + year).eq(1),
    year
  );
});
qs0AnomalyYear = qs0AnomalyYear
  .updateMask(quasiStable0)
  .rename('qs0_anomaly_year')
  .clip(roi);

// Paleta Viridis D — 9 colores para los 9 años (2017–2025)
var viridisPalette = [
  '#440154', '#472d7b', '#3b528b', '#2c728e', '#21908d',
  '#27ad81', '#5cc863', '#aadc32', '#fde725'
];

Map.addLayer(qs1AnomalyYear,
  {min: 2017, max: 2025, palette: viridisPalette},
  'QS1: año atípico sin retama (pred=0 único)', false);

Map.addLayer(qs0AnomalyYear,
  {min: 2017, max: 2025, palette: viridisPalette},
  'QS0: año atípico con retama (pred=1 único)', false);

// --- 6b. Inspector: serie temporal RF1 por pixel clickeado (consola) ---
Map.onClick(function(coords) {
  var point = ee.Geometry.Point([coords.lon, coords.lat]);

  predStack.reduceRegion({
    reducer : ee.Reducer.first(),
    geometry: point,
    scale   : 10
  }).evaluate(function(result) {
    var rows = exportYears.map(function(year) {
      var val = result['pred_' + year];
      return {
        year: year,
        RF1 : (val !== null && val !== undefined) ? val : 'NA'
      };
    });
    print('RF1 @ lon:' + coords.lon.toFixed(5) + ' lat:' + coords.lat.toFixed(5), rows);
  });
});


// =============================================================================
// 7. EXPORT: 4-category layer (1=s1, 2=qs1, 3=s0, 4=qs0)
// -----------------------------------------------------------------------------
// Exporting the categorized layer lets 04_RF2sampleExport.js decide at
// sampling time which categories to include, without re-running this script.
// Values: 1 = stable-retama (sumOnes==N)
//         2 = quasi-stable retama (sumOnes==N-1)
//         3 = stable-background (sumOnes==0)
//         4 = quasi-stable background (sumOnes==1)
// Asset name prefix: 03_
// =============================================================================
Export.image.toAsset({
  image      : stableLayer,
  description: 'Export_03_StableCategories',
  assetId    : ASSET_PREFIX + '03_RF2_stable_categories_rf1patchFiltered',
  region     : roi,
  scale      : 10,
  maxPixels  : 1e13
});

if (EXPORT_TO_DRIVE) {
  Export.image.toDrive({
    image         : stableLayer,
    description   : 'Drive_03_StableCategories',
    folder        : DRIVE_FOLDER,
    fileNamePrefix: '03_RF2_stable_categories_rf1patchFiltered',
    region        : roi,
    scale         : 10,
    maxPixels     : 1e13,
    fileFormat    : 'GeoTIFF'
  });
}
