// =============================================================================
// 07: RF2 — Patch filter on raw predictions → final RF2 assets
// -----------------------------------------------------------------------------
// Prerequisite : assets RF2_raw_prediction_YYYY  (from 06_RF2predict.js)
// Produces     : assets RF2_prediction_YYYY  (uint8, patch-filtered)
// -----------------------------------------------------------------------------
// Applying connectedPixelCount on already-materialised rasters is much cheaper
// than applying it inside the classification graph: GEE only needs to solve the
// neighbourhood topology on a simple binary image, not re-run the RF classifier
// for each tile's extended neighbourhood.
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX  = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER  = 'Retamap/Retamap_GEE_Exports';
var exportYears   = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var roi           = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');
var MIN_PATCH_PIX = 5;   // minimum connected-pixel count for class-1 patches
var MAX_NEIGHBORS = 50;  // neighbourhood radius for connectedPixelCount
var RUN_SUFFIX    = '_QS_n10k_compSamp';

// =============================================================================
// 1. LOAD RAW PREDICTIONS
// =============================================================================
var rawCollection = ee.ImageCollection.fromImages(
  exportYears.map(function(year) {
    return ee.Image(ASSET_PREFIX + 'RF2_raw_prediction_' + year + RUN_SUFFIX).set('year', year);
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

var filteredCollection = rawCollection.map(applyPatchFilter);

// =============================================================================
// 3. VISUALIZATION
// =============================================================================
Map.centerObject(roi, 10);

exportYears.forEach(function(year) {
  // Raw (before filter) — useful for comparison
  Map.addLayer(
    ee.Image(rawCollection.filter(ee.Filter.eq('year', year)).first()).clip(roi),
    {min: 0, max: 1, palette: ['#d3d3d3', '#2d9c00']},
    'RF2 raw ' + year, false
  );
  // Filtered (final)
  Map.addLayer(
    ee.Image(filteredCollection.filter(ee.Filter.eq('year', year)).first()).clip(roi),
    {min: 0, max: 1, palette: ['#d3d3d3', '#2d9c00']},
    'RF2 ' + year, false
  );
});

// =============================================================================
// 4. EXPORT FINAL PREDICTIONS PER YEAR
// =============================================================================
exportYears.forEach(function(year) {
  var img = ee.Image(
    filteredCollection.filter(ee.Filter.eq('year', year)).first()
  ).select('classification').clip(roi).uint8();

  Export.image.toAsset({
    image      : img,
    description: 'Export_RF2_' + year + RUN_SUFFIX,
    assetId    : ASSET_PREFIX + 'RF2_prediction_' + year + RUN_SUFFIX,
    region     : roi,
    scale      : 10,
    maxPixels  : 1e13
  });

  Export.image.toDrive({
    image          : img,
    description    : 'Drive_RF2_' + year + RUN_SUFFIX,
    folder         : DRIVE_FOLDER,
    fileNamePrefix : 'RF2_prediction_' + year + RUN_SUFFIX,
    region         : roi,
    scale          : 10,
    maxPixels      : 1e13,
    fileFormat     : 'GeoTIFF'
  });
});
