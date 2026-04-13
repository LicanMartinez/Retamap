// =============================================================================
// 04: RF2 — Stratified sampling from stable regions + Complementary samples
// =============================================================================
// Strategy:
//   1. bounds() on region: eliminates point-in-polygon over 14k vertices
//   2. Two-step: stratifiedSample over a single static-band image (once),
//      then sampleRegions per year (point lookups, not a spatial scan)
//   3. Step 1 at scale:10 over 1 band → pixel centres align with the native
//      Sentinel-2 grid, no ambiguity in spectral extraction
//   4. One samplingPoints task + 9 sampleRegions tasks (parallel, lightweight)
//   5. Complementary samples: sampleRegions over samplingLabel inside polygons
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX  = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER  = 'Retamap/Retamap_GEE_Exports';
var exportYears   = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var roi           = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');
var BANDS         = ['B2', 'B3', 'B4', 'B8', 'NDYI',
                     'B2_feb', 'B3_feb', 'B4_feb', 'B8_feb', 'NDYI_feb'];
var N_PER_CLASS   = 20000;
var RUN_SUFFIX    = '_noQS_n10k_compSamp';

var INCLUDE_S1  = true;
var INCLUDE_QS1 = false;
var INCLUDE_S0  = true;
var INCLUDE_QS0 = false;

var USE_COMPLEMENTARY      = true;
var COMPLEMENTARY_POLYGONS = ASSET_PREFIX + 'RF2_complementarySamplingStableRegions_kmsBari';

// =============================================================================
// 1. LOAD ASSETS
// =============================================================================
var stableCat = ee.Image(ASSET_PREFIX + 'RF2_stable_categories').select('stable_class');

var samplingRegion = roi.geometry().bounds();

// =============================================================================
// 2. BUILD SAMPLING LABEL
// =============================================================================
var mask1 = ee.Image(0).byte();
var mask0 = ee.Image(0).byte();

if (INCLUDE_S1)  mask1 = mask1.or(stableCat.eq(1));
if (INCLUDE_QS1) mask1 = mask1.or(stableCat.eq(2));
if (INCLUDE_S0)  mask0 = mask0.or(stableCat.eq(3));
if (INCLUDE_QS0) mask0 = mask0.or(stableCat.eq(4));

var samplingLabel = ee.Image(0).byte()
  .where(mask0, 0)
  .where(mask1, 1)
  .updateMask(mask1.or(mask0))
  .rename('stable_label');

// =============================================================================
// 3. STEP 1: SAMPLING LOCATIONS (1 band, scale:10, executed once)
// =============================================================================
var samplePoints = samplingLabel.stratifiedSample({
  numPoints   : N_PER_CLASS,
  classBand   : 'stable_label',
  region      : samplingRegion,
  scale       : 10,
  classValues : [0, 1],
  classPoints : [N_PER_CLASS, N_PER_CLASS],
  geometries  : true,
  tileScale   : 16
});

// =============================================================================
// 3.5. COMPLEMENTARY SAMPLES (sampleRegions within polygons)
// -----------------------------------------------------------------------------
// Extracts samplingLabel values only for unmasked pixels that fall inside
// the compPolys geometries. These points are merged into samplePoints before
// the per-year spectral extraction (Step 2), enriching coverage in areas
// that stratifiedSample may under-represent.
// =============================================================================
if (USE_COMPLEMENTARY) {
  var compPolys = ee.FeatureCollection(COMPLEMENTARY_POLYGONS);

  var compPoints = samplingLabel.sampleRegions({
    collection  : compPolys,
    scale       : 10,
    geometries  : true,
    tileScale   : 16
  });

  samplePoints = samplePoints.merge(compPoints);
}

// =============================================================================
// 4. STEP 2: SPECTRAL EXTRACTION PER YEAR (sampleRegions, parallel tasks)
// =============================================================================
exportYears.forEach(function(year) {
  var samples = ee.Image(ASSET_PREFIX + 'MergedBands_' + year)
    .select(BANDS)
    .sampleRegions({
      collection  : samplePoints,
      properties  : ['stable_label'],
      scale       : 10,
      geometries  : true,
      tileScale   : 16
    })
    .map(function(f) { return f.set('year', year); });

  Export.table.toAsset({
    collection  : samples,
    description : 'Export_RF2_samples_' + year + RUN_SUFFIX,
    assetId     : ASSET_PREFIX + 'RF2_samples_' + year + RUN_SUFFIX
  });

  Export.table.toDrive({
    collection     : samples,
    description    : 'Drive_RF2_samples_' + year + RUN_SUFFIX,
    folder         : DRIVE_FOLDER,
    fileNamePrefix : 'RF2_samples_' + year + RUN_SUFFIX,
    fileFormat     : 'CSV'
  });
});
