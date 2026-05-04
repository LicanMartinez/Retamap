// =============================================================================
// 08: RF2 — Patch filter on gap-filled predictions → final RF2 assets
// -----------------------------------------------------------------------------
// Prerequisite : assets 07_RF2_gapFill_YYYY{SUFFIX}  (from 07_RF2gapFill.js)
// Produces     : assets 08_RF2_prediction_YYYY{SUFFIX}  (uint8, patch-filtered)
// -----------------------------------------------------------------------------
// Applying connectedPixelCount on already-materialised rasters is much cheaper
// than applying it inside the classification graph: GEE only needs to solve the
// neighbourhood topology on a simple binary image, not re-run the RF classifier
// for each tile's extended neighbourhood.
// Input changed from 06_RF2_raw to 07_RF2_gapFill so the patch filter works
// on temporally-corrected predictions.
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER    = 'Retamap_GEE_Exports';
var EXPORT_TO_DRIVE = false;  // toggle: true → also export to Google Drive

var exportYears   = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var roi           = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');
var MIN_PATCH_PIX = 3;   // minimum connected-pixel count for class-1 patches
var MAX_NEIGHBORS = 50;  // neighbourhood radius for connectedPixelCount
var RUN_SUFFIX    = '_trimmedRF1comp_s1n12k_qs1n12k_s0n12k_qs0n12k';

// =============================================================================
// 1. LOAD GAP-FILLED PREDICTIONS (07_RF2_gapFill_YYYY)
// =============================================================================
var gapFillCollection = ee.ImageCollection.fromImages(
  exportYears.map(function(year) {
    return ee.Image(ASSET_PREFIX + '07_RF2_gapFill_' + year + RUN_SUFFIX).set('year', year);
  })
);

// =============================================================================
// 2. PATCH FILTER FUNCTION
// -----------------------------------------------------------------------------
// Class-1 patches with fewer than MIN_PATCH_PIX connected pixels → masked to
// nodata. Class-0 pixels are preserved unchanged.
// =============================================================================
var applyPatchFilter = function(image) {
  var pred = image.select('classification');

  var connected = pred
    .selfMask()                              // isolate class-1 pixels
    .connectedPixelCount(MAX_NEIGHBORS, true);

  // Keep class-1 only where patch is large enough; class-0 always passes through
  var class1_filtered = pred
    .updateMask(connected.gte(MIN_PATCH_PIX))
    .where(pred.eq(0), 0);                   // restore class-0 under the mask

  // Merge: class-0 from original, filtered class-1
  return pred
    .where(pred.eq(1), class1_filtered)
    .updateMask(
      pred.eq(0).or(connected.gte(MIN_PATCH_PIX))
    )
    .rename('classification')
    .set('year', image.get('year'));
};

var filteredCollection = gapFillCollection.map(applyPatchFilter);

// =============================================================================
// 3. VISUALIZATION — gap-filled vs final (patch-filtered)
// =============================================================================
Map.centerObject(roi, 10);

exportYears.forEach(function(year) {
  // Gap-filled (before patch filter) — useful for comparison
  Map.addLayer(
    ee.Image(gapFillCollection.filter(ee.Filter.eq('year', year)).first()).clip(roi),
    {min: 0, max: 1, palette: ['#d3d3d3', '#2d9c00']},
    'RF2 gap-fill ' + year, false
  );
  // Final (patch-filtered)
  Map.addLayer(
    ee.Image(filteredCollection.filter(ee.Filter.eq('year', year)).first()).clip(roi),
    {min: 0, max: 1, palette: ['#d3d3d3', '#CC0000']},
    'RF2 final ' + year, false
  );
});

// =============================================================================
// 4. EXPORT FINAL PREDICTIONS PER YEAR (08_RF2_prediction_YYYY{SUFFIX})
// =============================================================================
exportYears.forEach(function(year) {
  var img = ee.Image(
    filteredCollection.filter(ee.Filter.eq('year', year)).first()
  ).select('classification').clip(roi).uint8();

  Export.image.toAsset({
    image      : img,
    description: 'Export_08_RF2_' + year + RUN_SUFFIX,
    assetId    : ASSET_PREFIX + '08_RF2_prediction_' + year + RUN_SUFFIX,
    region     : roi,
    scale      : 10,
    maxPixels  : 1e13
  });

  if (EXPORT_TO_DRIVE) {
    Export.image.toDrive({
      image         : img,
      description   : 'Drive_08_RF2_' + year + RUN_SUFFIX,
      folder        : DRIVE_FOLDER,
      fileNamePrefix: '08_RF2_prediction_' + year + RUN_SUFFIX,
      region        : roi,
      scale         : 10,
      maxPixels     : 1e13,
      fileFormat    : 'GeoTIFF'
    });
  }
});
