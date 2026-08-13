// =============================================================================
// 08: RF2 — Patch filter on gap-filled predictions → final RF2 assets
// -----------------------------------------------------------------------------
// Role         : MAP PRODUCTION, step 8 of 8 (final map). Section 2.2.5.
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
var RUN_SUFFIX    = '_trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k';

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
// nodata.
//
// ⚠ The output is a PRESENCE-ONLY asset: class-0 pixels end up masked too,
// despite what the .where(pred.eq(0), 0) below suggests. `connected` is derived
// from pred.selfMask(), so it is masked wherever pred == 0; EE's .or() does not
// short-circuit, it INTERSECTS masks, so the updateMask() at the end keeps only
// "class 1 in a patch of ≥ MIN_PATCH_PIX". Everything else is nodata, NOT 0.
//
// This is load-bearing downstream and must not be "fixed" casually:
//   09v_valGenerator.js  §2  unmask(0)s the maps and takes the valid-territory
//                            mask from 01_MergedBands instead of from here.
//   09v_valMapLabels.js      leaves the labels masked on purpose (blank → 0 in R).
// It also means any area statistic has to reintroduce the observation footprint
// from 01_MergedBands — see 13_areaAudit.js.
// =============================================================================
var applyPatchFilter = function(image) {
  var pred = image.select('classification');

  var connected = pred
    .selfMask()                              // isolate class-1 pixels
    .connectedPixelCount(MAX_NEIGHBORS, true);

  // Keep class-1 only where the patch is large enough.
  var class1_filtered = pred
    .updateMask(connected.gte(MIN_PATCH_PIX))
    .where(pred.eq(0), 0);

  // NOTE: the mask below intersects with `connected`'s mask (see header), so the
  // result is presence-only — class-0 does NOT survive this step.
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
