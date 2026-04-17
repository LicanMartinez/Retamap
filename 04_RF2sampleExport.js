// =============================================================================
// 04: RF2 — Stratified sampling from stable regions + Complementary samples
//          + Train/Val split for holdout validation
// =============================================================================
// Strategy:
//   1. bounds() on region: eliminates point-in-polygon over 14k vertices
//   2. Two-step: stratifiedSample over a single static-band image (once),
//      then sampleRegions per year (point lookups, not a spatial scan)
//   3. Step 1 at scale:10 over 1 band → pixel centres align with the native
//      Sentinel-2 grid, no ambiguity in spectral extraction
//   4. One stratifiedSample task + 18 sampleRegions tasks (9 _train + 9 _val)
//   5. Complementary samples: sampleRegions over binaryLabel inside polygons
//   6. MOD 1: per-subclass sample counts (N_S1/N_QS1/N_S0/N_QS0)
//   7. MOD 5: random train/val split before spectral extraction
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER    = 'Retamap/Retamap_GEE_Exports';
var EXPORT_TO_DRIVE = false;  // toggle: true → also export to Google Drive

var exportYears = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var roi         = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');
var BANDS       = ['B2', 'B3', 'B4', 'B8', 'NDYI',
                   'B2_feb', 'B3_feb', 'B4_feb', 'B8_feb', 'NDYI_feb'];

// ── Subclass sample counts (MOD 1) ──────────────────────────────────────────
// Set a count to 0 to exclude that subclass entirely.
var N_S1  = 15000;  // stable retama       (category 1)
var N_QS1 = 15000;      // quasi-stable retama (category 2) — 0 to exclude
var N_S0  = 15000;  // stable background   (category 3)
var N_QS0 = 15000;      // quasi-stable bg     (category 4) — 0 to exclude

// Suffix template — must reflect the counts above.
var RUN_SUFFIX = '_s1n15k_qs1n15k_s0n15k_qs0n15k';

// ── Train / validation split (MOD 5) ────────────────────────────────────────
var VAL_FRACTION = 0.2;   // fraction of sample points held out for validation
var RANDOM_SEED  = 42;

// ── Complementary samples ───────────────────────────────────────────────────
var USE_COMPLEMENTARY      = false;
var COMPLEMENTARY_POLYGONS = ASSET_PREFIX + 'RF2_complementarySamplingStableRegions_kmsBari';

// =============================================================================
// 1. LOAD ASSETS
// =============================================================================
var stableCat = ee.Image(ASSET_PREFIX + '03_RF2_stable_categories').select('stable_class');

var samplingRegion = roi.geometry().bounds();

// =============================================================================
// 2. BUILD classValues / classPoints (skipping categories with count == 0)
// =============================================================================
var classValues = [];
var classPoints = [];
if (N_S1  > 0) { classValues.push(1); classPoints.push(N_S1);  }
if (N_QS1 > 0) { classValues.push(2); classPoints.push(N_QS1); }
if (N_S0  > 0) { classValues.push(3); classPoints.push(N_S0);  }
if (N_QS0 > 0) { classValues.push(4); classPoints.push(N_QS0); }

// =============================================================================
// 3. STEP 1: SAMPLING LOCATIONS (1 band, scale:10, executed once)
// =============================================================================
var samplePoints = stableCat.stratifiedSample({
  numPoints  : 0,
  classBand  : 'stable_class',
  region     : samplingRegion,
  scale      : 10,
  classValues: classValues,
  classPoints: classPoints,
  geometries : true,
  tileScale  : 16
});

// Remap 4 categories → binary label for downstream RF2 training compatibility:
// categories 1 & 2 (retama) → 1,  categories 3 & 4 (background) → 0
samplePoints = samplePoints.map(function(f) {
  var sc = ee.Number(f.get('stable_class'));
  return f.set('stable_label', ee.Algorithms.If(sc.lte(2), 1, 0));
});

// =============================================================================
// 3.5. COMPLEMENTARY SAMPLES (sampleRegions within polygons)
// -----------------------------------------------------------------------------
// Uses a binary label derived from stableCat for compatibility with Step 2.
// =============================================================================
if (USE_COMPLEMENTARY) {
  var compPolys = ee.FeatureCollection(COMPLEMENTARY_POLYGONS);

  // Binary label for complementary sampleRegions
  var binaryLabel = stableCat.remap([1, 2, 3, 4], [1, 1, 0, 0]).rename('stable_label');

  var compPoints = binaryLabel.sampleRegions({
    collection: compPolys,
    scale     : 10,
    geometries: true,
    tileScale : 16
  });

  samplePoints = samplePoints.merge(compPoints);
}


// =============================================================================
// 4. TRAIN / VALIDATION SPLIT (MOD 5)
// -----------------------------------------------------------------------------
// Split happens BEFORE spectral extraction so the same spatial points are used
// consistently across all 9 years. Each year produces two exports: _train and
// _val. Script 05 uses only _train; script 09 uses _val for holdout evaluation.
// =============================================================================
samplePoints = samplePoints.randomColumn('rand', RANDOM_SEED);
var trainPoints = samplePoints.filter(ee.Filter.gte('rand', VAL_FRACTION));
var valPoints   = samplePoints.filter(ee.Filter.lt('rand', VAL_FRACTION));

// =============================================================================
// 4.5. MAP VISUALIZATION OF SAMPLE POINTS
// =============================================================================
Map.centerObject(roi, 11);

// Visualize by the 4 specific categories
Map.addLayer(samplePoints.filter(ee.Filter.eq('stable_class', 1)), {color: '#d62728'}, 'Class 1: Stable Retama');
Map.addLayer(samplePoints.filter(ee.Filter.eq('stable_class', 2)), {color: '#ff7f0e'}, 'Class 2: Quasi-stable Retama');
Map.addLayer(samplePoints.filter(ee.Filter.eq('stable_class', 3)), {color: '#1f77b4'}, 'Class 3: Stable Background');
Map.addLayer(samplePoints.filter(ee.Filter.eq('stable_class', 4)), {color: '#17becf'}, 'Class 4: Quasi-stable Background');

// Visualize by the binary label (includes complementary samples). Hidden by default.
Map.addLayer(samplePoints.filter(ee.Filter.eq('stable_label', 1)), {color: 'red'}, 'Binary Label: Retama (1)', false);
Map.addLayer(samplePoints.filter(ee.Filter.eq('stable_label', 0)), {color: 'blue'}, 'Binary Label: Background (0)', false);



// =============================================================================
// 5. STEP 2: SPECTRAL EXTRACTION PER YEAR (sampleRegions, parallel tasks)
// =============================================================================
exportYears.forEach(function(year) {
  var img = ee.Image(ASSET_PREFIX + '01_MergedBands_' + year).select(BANDS);

  // ── Train samples ──────────────────────────────────────────────────────────
  var trainSamples = img.sampleRegions({
    collection: trainPoints,
    properties: ['stable_label'],
    scale     : 10,
    geometries: true,
    tileScale : 16
  }).map(function(f) { return f.set('year', year); });

  Export.table.toAsset({
    collection : trainSamples,
    description: 'Export_04_RF2_samples_' + year + RUN_SUFFIX + '_train',
    assetId    : ASSET_PREFIX + '04_RF2_samples_' + year + RUN_SUFFIX + '_train'
  });

  if (EXPORT_TO_DRIVE) {
    Export.table.toDrive({
      collection    : trainSamples,
      description   : 'Drive_04_RF2_samples_' + year + RUN_SUFFIX + '_train',
      folder        : DRIVE_FOLDER,
      fileNamePrefix: '04_RF2_samples_' + year + RUN_SUFFIX + '_train',
      fileFormat    : 'CSV'
    });
  }

  // ── Val samples ────────────────────────────────────────────────────────────
  var valSamples = img.sampleRegions({
    collection: valPoints,
    properties: ['stable_label'],
    scale     : 10,
    geometries: true,
    tileScale : 16
  }).map(function(f) { return f.set('year', year); });

  Export.table.toAsset({
    collection : valSamples,
    description: 'Export_04_RF2_samples_' + year + RUN_SUFFIX + '_val',
    assetId    : ASSET_PREFIX + '04_RF2_samples_' + year + RUN_SUFFIX + '_val'
  });

  if (EXPORT_TO_DRIVE) {
    Export.table.toDrive({
      collection    : valSamples,
      description   : 'Drive_04_RF2_samples_' + year + RUN_SUFFIX + '_val',
      folder        : DRIVE_FOLDER,
      fileNamePrefix: '04_RF2_samples_' + year + RUN_SUFFIX + '_val',
      fileFormat    : 'CSV'
    });
  }
});
