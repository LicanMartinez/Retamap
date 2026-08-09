// =============================================================================
// SIDE EXPLORATION — NOT PART OF THE CANONICAL PIPELINE.
// Explored, measured against the 09v validation, and NOT ADOPTED (2026-08-09).
// Full write-up and numbers in ./README.md. Canonical: 08_RF2patchFilter.js in
// ../../. NOTE the globalRF2_noGapFill arm below is still of independent
// interest: comparing it against the published map isolates what gap-fill
// contributes, and the 09v comparison found that to be very little.
// =============================================================================
// 08y: RF2 VARIANT (yearlyRF2) — patch filter on the variant AND on a MATCHED
//      no-gap-fill baseline of the canonical global RF2
// -----------------------------------------------------------------------------
// Prerequisite : assets 06y_RF2_raw_prediction_YYYY_yearlyRF2      (from 06y)
//              + assets 06_RF2_raw_prediction_YYYY{CANON_SUFFIX}   (ALREADY
//                exported by the canonical run — nothing to re-run)
// Produces     : 08y_RF2_prediction_YYYY_yearlyRF2_noGapFill
//                08y_RF2_prediction_YYYY_globalRF2_noGapFill
// -----------------------------------------------------------------------------
// WHY THE BASELINE EXISTS
// Script 07 (gap-fill) is a TEMPORAL filter: `sumOnes == N-1 → 1` and
// `sumOnes == 1 → 0`. With N == 2 both thresholds collapse onto 1, so every
// pixel where 2017 and 2025 disagree gets SWAPPED — the change signal is
// destroyed. Gap-fill has no meaningful 2-year analogue, so the variant skips
// it entirely.
//
// That means the variant cannot be compared against the PUBLISHED canonical map
// without conflating two changes (per-year RF2 *and* no gap-fill) — and since
// gap-fill mainly repairs omissions, exactly where the canonical map is weakest,
// that confound would bias against the variant. So this script also runs the
// SAME patch-filter function over the canonical raw predictions, producing a
// matched `globalRF2_noGapFill` baseline. The A/B is then identical in every
// respect except RF2 training.
//
// IF THIS IS EVER RUN FOR ALL 9 YEARS
// Add a 07y_RF2gapFill.js (calque of 07 with the 06y_/07y_ stems and the full
// exportYears), then flip the yearlyRF2 entry of RUNS below to
// inStem: '07y_RF2_gapFill_' and outSuffix: '_yearlyRF2' (no marker, matching
// the canonical convention). One line of config.
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER    = 'Retamap/Retamap_GEE_Exports';
var EXPORT_TO_DRIVE = false;  // toggle: true → also export to Google Drive

var exportYears   = [2017, 2025];
var roi           = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');

// Identical to canonical 08 — the filter must not be a difference between arms.
var MIN_PATCH_PIX = 3;   // minimum connected-pixel count for class-1 patches
var MAX_NEIGHBORS = 50;  // neighbourhood radius for connectedPixelCount

var CANON_SUFFIX = '_trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k';

// The two arms of the experiment. Both go through applyPatchFilter() below, so
// the filtering is guaranteed to be the same operation, not a re-implementation.
var RUNS = [
  { label    : 'yearly RF2',
    inStem   : '06y_RF2_raw_prediction_',
    inSuffix : '_yearlyRF2',
    outSuffix: '_yearlyRF2_noGapFill',
    palette  : '#1f6fd0' },
  { label    : 'global RF2',
    inStem   : '06_RF2_raw_prediction_',
    inSuffix : CANON_SUFFIX,
    outSuffix: '_globalRF2_noGapFill',
    palette  : '#e08214' }
];

// =============================================================================
// 1. PATCH FILTER FUNCTION  (verbatim from 08_RF2patchFilter.js)
// -----------------------------------------------------------------------------
// Class-1 patches with fewer than MIN_PATCH_PIX connected pixels → masked to
// nodata. Class-0 pixels are preserved unchanged. Year-count agnostic: it maps
// over each image independently.
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

// =============================================================================
// 2. BUILD FILTERED IMAGES FOR BOTH ARMS
// =============================================================================
var filtered = {};   // filtered[outSuffix][year] = image

RUNS.forEach(function(run) {
  filtered[run.outSuffix] = {};
  exportYears.forEach(function(year) {
    var raw = ee.Image(ASSET_PREFIX + run.inStem + year + run.inSuffix)
      .set('year', year);
    filtered[run.outSuffix][year] = applyPatchFilter(raw);
  });
});

// =============================================================================
// 3. VISUALIZATION — published canonical vs the two matched no-gap-fill arms
// -----------------------------------------------------------------------------
// All off by default. Turning the three on together shows where the arms differ:
// canonical-vs-globalNoGF isolates the gap-fill effect, globalNoGF-vs-yearlyNoGF
// isolates the per-year training effect.
// =============================================================================
Map.centerObject(roi, 10);

exportYears.forEach(function(year) {
  // Published canonical map (WITH gap-fill) — context only, not produced here
  Map.addLayer(
    ee.Image(ASSET_PREFIX + '08_RF2_prediction_' + year + CANON_SUFFIX)
      .select('classification').selfMask().clip(roi),
    {min: 0, max: 1, palette: ['#CC0000']},
    'canonical + gapFill ' + year, false
  );

  RUNS.forEach(function(run) {
    Map.addLayer(
      filtered[run.outSuffix][year].selfMask().clip(roi),
      {min: 0, max: 1, palette: [run.palette]},
      run.label + ' noGapFill ' + year, false
    );
  });
});

// =============================================================================
// 4. EXPORT FINAL PREDICTIONS (08y_RF2_prediction_YYYY{outSuffix})
// =============================================================================
RUNS.forEach(function(run) {
  exportYears.forEach(function(year) {
    var img = filtered[run.outSuffix][year]
      .select('classification').clip(roi).uint8();

    Export.image.toAsset({
      image      : img,
      description: 'Export_08y_RF2_' + year + run.outSuffix,
      assetId    : ASSET_PREFIX + '08y_RF2_prediction_' + year + run.outSuffix,
      region     : roi,
      scale      : 10,
      maxPixels  : 1e13
    });

    if (EXPORT_TO_DRIVE) {
      Export.image.toDrive({
        image         : img,
        description   : 'Drive_08y_RF2_' + year + run.outSuffix,
        folder        : DRIVE_FOLDER,
        fileNamePrefix: '08y_RF2_prediction_' + year + run.outSuffix,
        region        : roi,
        scale         : 10,
        maxPixels     : 1e13,
        fileFormat    : 'GeoTIFF'
      });
    }
  });
});
