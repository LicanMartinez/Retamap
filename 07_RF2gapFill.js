// =============================================================================
// 07: RF2 — Temporal gap-fill of raw predictions
// -----------------------------------------------------------------------------
// Prerequisite : assets 06_RF2_raw_prediction_YYYY{SUFFIX}  (from 06_RF2predict.js)
// Produces     : assets 07_RF2_gapFill_YYYY{SUFFIX}  (uint8)
// -----------------------------------------------------------------------------
// Logic (MOD 2):
//   Stack all raw RF2 predictions and compute sumOnes per pixel.
//   For each year:
//     - quasi-stable retama  (sumOnes == N-1 AND this year pred==0) → fill to 1
//     - quasi-stable bg      (sumOnes == 1   AND this year pred==1) → fill to 0
//   All other pixels pass through unchanged.
//
//   This corrects one-year anomalies caused by cloud / atmospheric noise
//   before the patch filter (08_RF2patchFilter.js).
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER    = 'Retamap/Retamap_GEE_Exports';
var EXPORT_TO_DRIVE = false;  // toggle: true → also export to Google Drive

var exportYears = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var roi         = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');
var RUN_SUFFIX  = '_s1.10k_qs1.0_s0.20k_qs0.0';

// =============================================================================
// 1. LOAD RAW PREDICTIONS AND BUILD MULTI-BAND STACK
// =============================================================================
var predList = exportYears.map(function(year) {
  return ee.Image(ASSET_PREFIX + '06_RF2_raw_prediction_' + year + RUN_SUFFIX)
    .select('classification')
    .rename('pred_' + year);
});

var predStack = ee.Image(predList);
var N = exportYears.length;  // 9 years

// Per-pixel count of years with RF2 raw == 1
var sumOnes = predStack.reduce(ee.Reducer.sum()).rename('sum_ones');

// =============================================================================
// 2. DIAGNOSTIC: pixel counts for gap-fill targets
// =============================================================================
var diagMask_qs1 = sumOnes.eq(N - 1);  // pixels that are quasi-stable retama
var diagMask_qs0 = sumOnes.eq(1);      // pixels that are quasi-stable background

print('Gap-fill diagnostic:');
print('  Quasi-stable retama pixels (sumOnes == ' + (N - 1) + '):',
  diagMask_qs1.reduceRegion({
    reducer: ee.Reducer.sum(), geometry: roi.geometry(),
    scale: 10, maxPixels: 1e13, tileScale: 16
  }).get('sum_ones_sum')
);
print('  Quasi-stable bg pixels (sumOnes == 1):',
  diagMask_qs0.reduceRegion({
    reducer: ee.Reducer.sum(), geometry: roi.geometry(),
    scale: 10, maxPixels: 1e13, tileScale: 16
  }).get('sum_ones_sum')
);

// =============================================================================
// 3. GAP-FILL PER YEAR AND EXPORT
// =============================================================================
exportYears.forEach(function(year) {
  var original = predStack.select('pred_' + year);

  var gapFilled = original
    // quasi-stable retama: N-1 years == 1, this year == 0 → fill with 1
    .where(sumOnes.eq(N - 1).and(original.eq(0)), 1)
    // quasi-stable background: 1 year == 1, this year == 1 → fill with 0
    .where(sumOnes.eq(1).and(original.eq(1)), 0)
    .rename('classification')
    .clip(roi)
    .uint8();

  Export.image.toAsset({
    image      : gapFilled,
    description: 'Export_07_RF2_gapFill_' + year + RUN_SUFFIX,
    assetId    : ASSET_PREFIX + '07_RF2_gapFill_' + year + RUN_SUFFIX,
    region     : roi,
    scale      : 10,
    maxPixels  : 1e13
  });

  if (EXPORT_TO_DRIVE) {
    Export.image.toDrive({
      image         : gapFilled,
      description   : 'Drive_07_RF2_gapFill_' + year + RUN_SUFFIX,
      folder        : DRIVE_FOLDER,
      fileNamePrefix: '07_RF2_gapFill_' + year + RUN_SUFFIX,
      region        : roi,
      scale         : 10,
      maxPixels     : 1e13,
      fileFormat    : 'GeoTIFF'
    });
  }
});

// =============================================================================
// 4. VISUALIZATION — raw vs gap-filled comparison
// =============================================================================
Map.centerObject(roi, 10);

exportYears.forEach(function(year) {
  Map.addLayer(
    predStack.select('pred_' + year).selfMask().clip(roi),
    {min: 1, max: 1, palette: ['#FF8C00']},
    'RF2 raw ' + year, false
  );
});

Map.addLayer(
  sumOnes.clip(roi),
  {min: 0, max: N, palette: ['white', 'yellow', 'orange', 'red', 'darkred']},
  'sumOnes (RF2 raw)', false
);
