// =============================================================================
// SIDE EXPLORATION — NOT PART OF THE CANONICAL PIPELINE.
// Explored, measured against the 09v validation, and NOT ADOPTED (2026-08-09):
// per-year training worsens retama omission in both validated years and breaks
// the map's zero-commission property in the far background. Full write-up and
// numbers in ./README.md. The canonical pipeline is 05_RF2train.js in ../../.
// =============================================================================
// 05y: RF2 VARIANT (yearlyRF2) — one classifier PER YEAR → classifier assets
// -----------------------------------------------------------------------------
// Prerequisite : assets 04_RF2_samples_YYYY{SAMPLES_SUFFIX}_train
//                (from 04_RF2sampleExport.js — reused unchanged, NOT re-exported)
// Produces     : assets 05y_RF2_classifier_YYYY{RUN_SUFFIX}
// -----------------------------------------------------------------------------
// EXPLORATORY VARIANT — the canonical pipeline (05_RF2train.js) is untouched.
//
// Canonical 05 pools the 9 per-year _train assets (~360k features) into ONE
// global RF2. This variant instead trains ONE RF2 per year on that year's
// samples only: the sampling LOCATIONS are identical across years (04 runs
// stratifiedSample once and then sampleRegions per year), so the only thing
// that changes between the per-year models is the spectral signal of that year.
// Hypothesis: per-year phenology / atmospheric conditions / mosaic quality make
// a year-specific model fit better than one global model.
//
// Run only for the two validated years (2017, 2025) — extend trainYears to the
// full 9 if the validation comparison says the variant is worth it.
//
// TWO SUFFIXES, on purpose:
//   SAMPLES_SUFFIX → reads the EXISTING 04_ sample assets (long canonical name)
//   RUN_SUFFIX     → writes the NEW variant assets (short readable name)
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX = 'projects/ee-licanemartinez/assets/Retamap/';

// Suffix of the sample assets being READ (canonical run — do not shorten, it is
// baked into the already-exported 04_ asset names).
var SAMPLES_SUFFIX = '_trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k';
// Suffix of the variant assets being WRITTEN.
var RUN_SUFFIX     = '_yearlyRF2';

var trainYears = [2017, 2025];   // validated years only; extend to the 9 later

// Must stay byte-identical to the canonical run: 10 bands, order matters for
// inputProperties. Together with numberOfTrees=500 these are the two constants
// that keep the variant comparable to the published map.
var BANDS      = ['B2', 'B3', 'B4', 'B8', 'NDYI', 'B2_feb', 'B3_feb', 'B4_feb', 'B8_feb', 'NDYI_feb'];
var N_TREES    = 500;

// =============================================================================
// 1. TRAIN ONE RF2 PER YEAR + EXPORT EACH CLASSIFIER
// -----------------------------------------------------------------------------
// Training on _train only (80%), same as canonical 05 — the _val split is a
// leftover of the removed internal holdout and is deliberately not used.
//
// Exporting the classifier (instead of training inline, as
// experimental/05b_RF2trainPerYear.js does) means 06y pays only the cost of
// classification per export task, and the per-year models stay inspectable.
// =============================================================================
trainYears.forEach(function(year) {
  var samples = ee.FeatureCollection(
    ASSET_PREFIX + '04_RF2_samples_' + year + SAMPLES_SUFFIX + '_train');

  print(year + ' — train samples:', samples.size());
  print(year + ' — class counts (stable_label):',
        samples.aggregate_histogram('stable_label'));

  var rf2 = ee.Classifier.smileRandomForest({numberOfTrees: N_TREES}).train({
    features       : samples,
    classProperty  : 'stable_label',
    inputProperties: BANDS
  });

  // Diagnostic for the experiment: does band importance shift between years?
  print(year + ' — variableImportance:',
        ee.Dictionary(rf2.explain().get('importance')));

  // ee.Classifier can only be exported to Asset (no Drive equivalent).
  Export.classifier.toAsset({
    classifier : rf2,
    description: 'Export_05y_RF2_classifier_' + year + RUN_SUFFIX,
    assetId    : ASSET_PREFIX + '05y_RF2_classifier_' + year + RUN_SUFFIX
  });
});
