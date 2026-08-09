// =============================================================================
// 06y: RF2 VARIANT (yearlyRF2) — load per-year classifier → raw predictions
// -----------------------------------------------------------------------------
// Prerequisite : assets 05y_RF2_classifier_YYYY{RUN_SUFFIX}  (from 05y)
//              + assets 01_MergedBands_YYYY                  (from 01)
// Produces     : assets 06y_RF2_raw_prediction_YYYY{RUN_SUFFIX}  (uint8)
// -----------------------------------------------------------------------------
// EXPLORATORY VARIANT — canonical 06_RF2predict.js is untouched.
//
// Each year is predicted with ITS OWN classifier and ITS OWN mosaic: the whole
// point of the variant. Structure is otherwise a calque of canonical 06.
//
// Next step is 08y_RF2patchFilter.js, NOT a gap-fill: script 07's temporal
// rules (sumOnes == N-1 → 1, sumOnes == 1 → 0) collapse onto each other when
// N == 2, which would swap every pixel where 2017 and 2025 disagree. Gap-fill
// has no valid 2-year analogue and is deliberately skipped here — see the note
// in 08y about what changes if this is ever run for all 9 years.
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER    = 'Retamap/Retamap_GEE_Exports';
var EXPORT_TO_DRIVE = false;  // toggle: true → also export to Google Drive

var RUN_SUFFIX  = '_yearlyRF2';
var exportYears = [2017, 2025];   // must match trainYears in 05y
var roi         = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');
var BANDS       = ['B2', 'B3', 'B4', 'B8', 'NDYI', 'B2_feb', 'B3_feb', 'B4_feb', 'B8_feb', 'NDYI_feb'];

// =============================================================================
// 1. PREDICT PER YEAR — one classifier per year
// -----------------------------------------------------------------------------
// .select(BANDS) is REQUIRED: 01_MergedBands_YYYY carries 11 bands (it also
// holds a `year` band from 01's addDateBands), and without the select the
// classifier would receive 11 inputs against 10 trained properties.
// =============================================================================
var rawByYear = {};

exportYears.forEach(function(year) {
  var rf2 = ee.Classifier.load(
    ASSET_PREFIX + '05y_RF2_classifier_' + year + RUN_SUFFIX);

  var img = ee.Image(ASSET_PREFIX + '01_MergedBands_' + year).select(BANDS);

  rawByYear[year] = img.classify(rf2).rename('classification').set('year', year);
});

// =============================================================================
// 2. VISUALIZATION
// =============================================================================
Map.centerObject(roi, 10);

exportYears.forEach(function(year) {
  Map.addLayer(
    rawByYear[year].clip(roi),
    {min: 0, max: 1, palette: ['#d3d3d3', '#2d9c00']},
    'RF2 yearly raw ' + year, false
  );
});

// =============================================================================
// 3. EXPORT RAW PREDICTIONS PER YEAR (06y_RF2_raw_prediction_YYYY{RUN_SUFFIX})
// =============================================================================
exportYears.forEach(function(year) {
  var img = rawByYear[year].select('classification').clip(roi).uint8();

  Export.image.toAsset({
    image      : img,
    description: 'Export_06y_RF2raw_' + year + RUN_SUFFIX,
    assetId    : ASSET_PREFIX + '06y_RF2_raw_prediction_' + year + RUN_SUFFIX,
    region     : roi,
    scale      : 10,
    maxPixels  : 1e13
  });

  if (EXPORT_TO_DRIVE) {
    Export.image.toDrive({
      image         : img,
      description   : 'Drive_06y_RF2raw_' + year + RUN_SUFFIX,
      folder        : DRIVE_FOLDER,
      fileNamePrefix: '06y_RF2_raw_prediction_' + year + RUN_SUFFIX,
      region        : roi,
      scale         : 10,
      maxPixels     : 1e13,
      fileFormat    : 'GeoTIFF'
    });
  }
});
