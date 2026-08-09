// =============================================================================
// 09v_valMapLabels: extract MAP labels for several map versions at the EXISTING
//                   validation points → Drive CSV for the R comparison
// -----------------------------------------------------------------------------
// Prerequisite : asset 09v_valPoints_fc  (blind FC, 919 features, only `pxid`)
//              + the 08/08y prediction assets listed in MAPS below
// Produces     : Drive CSV 09v_valMapLabels.csv  (pxid + one column per map/year)
// -----------------------------------------------------------------------------
// WHY THIS SCRIPT EXISTS — AND WHY NOT RE-RUN 09v_valGenerator.js
// The generator's stratifiedSample would draw a DIFFERENT set of points, which
// would orphan the ~1100 human labels already collected in the three Sheets.
// The point set, the strata, the area weights Wi and the reference labels are
// all frozen: the ONLY thing that changes when evaluating an alternative map is
// the map label at each point. That is exactly what this script recomputes.
//
// Evaluating a different map on this sample stays unbiased — strata and Wi are
// properties of the sampling DESIGN, not of the map — but power is uneven: any
// retama a variant adds OUTSIDE the canonical S1 is only sampled through S2
// (400 pts) and S3 (119 pts covering 93.5% of the ROI). The comparison report
// states this as a limitation.
//
// CONTROL COLUMN (the whole point of including the canonical map here)
// m<year>_canon is re-extracted from the very asset that produced the frozen
// m2017/m2025 columns of 09v_master_key.csv. The R checks the two agree row by
// row. If they do, the extraction is faithful and the variant columns can be
// trusted; if they don't, there is a grid misalignment and no metric is
// interpretable until it is resolved.
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER = 'Retamap/Retamap_GEE_Exports';

var FC_ASSET     = ASSET_PREFIX + '09v_valPoints_fc';
var VAL_YEARS    = [2017, 2025];
var CANON_SUFFIX = '_trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k';

// One entry per map version to evaluate. `tag` becomes the column suffix, so
// the R config (09v_valCompare.R) reads columns m<year>_<tag>. Evaluating a NEW
// map version = add an entry here and a matching one in CMP_VARIANTS.
//
// The two 08y_ entries below come from the `yearlyRF2` side exploration
// (explorations/yearlyRF2/, run and NOT adopted on 2026-08-09). They are kept
// because their assets still exist and the pair is a useful reference point:
// `canon` vs `globalNoGF` isolates what gap-fill contributes, independently of
// the per-year question. Drop or replace them freely.
var MAPS = [
  { tag: 'canon',
    stem: '08_RF2_prediction_',  suffix: CANON_SUFFIX },
  { tag: 'globalNoGF',
    stem: '08y_RF2_prediction_', suffix: '_globalRF2_noGapFill' },
  { tag: 'yearlyNoGF',
    stem: '08y_RF2_prediction_', suffix: '_yearlyRF2_noGapFill' }
];

// =============================================================================
// 1. POINTS — centroid of each blind 10 m square
// -----------------------------------------------------------------------------
// 09v_valRandomize.R built the FC polygons as 10 m squares in UTM19S centred on
// the sampled pixel centre, so the centroid recovers exactly the location the
// original m2017/m2025 were sampled at.
// =============================================================================
var pts = ee.FeatureCollection(FC_ASSET).map(function(f) {
  return ee.Feature(f.geometry().centroid(1), {pxid: f.get('pxid')});
});

print('Validation points loaded:', pts.size(), '(expected 919)');

// =============================================================================
// 2. BUILD THE MULTI-BAND LABEL STACK
// -----------------------------------------------------------------------------
// One band per (map version, year), named m<year>_<tag>. Left MASKED on purpose
// (no unmask): a pixel dropped by the patch filter has no value, and the R's
// zero_blank() maps the resulting blank cell to 0 — the correct semantics, and
// the same treatment the frozen master-key columns already received.
// =============================================================================
var bands = [];
MAPS.forEach(function(m) {
  VAL_YEARS.forEach(function(year) {
    bands.push(
      ee.Image(ASSET_PREFIX + m.stem + year + m.suffix)
        .select('classification')
        .rename('m' + year + '_' + m.tag)
    );
  });
});
var labelStack = ee.Image.cat(bands);

// Grid-A: the 4326 grid the maps were classified on (01_MergedBands exported at
// scale 10 with no crs). Reducing on this grid samples the very pixel the map
// assigned — same snap the validation App and the KMLs use.
var projA = ee.Image(ASSET_PREFIX + '01_MergedBands_' + VAL_YEARS[0])
              .select('B4').projection();

// =============================================================================
// 3. SAMPLE AT THE POINTS
// =============================================================================
var sampled = labelStack.reduceRegions({
  collection: pts,
  reducer   : ee.Reducer.first(),
  scale     : 10,
  crs       : projA,
  tileScale : 4
});

// Column order of the CSV: pxid, then m<year>_<tag> for every map/year.
var SELECTORS = ['pxid'];
MAPS.forEach(function(m) {
  VAL_YEARS.forEach(function(year) {
    SELECTORS.push('m' + year + '_' + m.tag);
  });
});
print('CSV columns:', SELECTORS);

// Quick console cross-check before the export finishes: per-map retama counts
// over the 919 points. canon vs globalNoGF should differ only slightly (the
// gap-fill effect); globalNoGF vs yearlyNoGF is the experiment.
SELECTORS.slice(1).forEach(function(col) {
  print('n(map=1) ' + col + ':',
        sampled.filter(ee.Filter.eq(col, 1)).size());
});

// =============================================================================
// 4. EXPORT → Drive CSV
// =============================================================================
Export.table.toDrive({
  collection    : sampled,
  description   : '09v_valMapLabels',
  folder        : DRIVE_FOLDER,
  fileNamePrefix: '09v_valMapLabels',
  fileFormat    : 'CSV',
  selectors     : SELECTORS
});

// =============================================================================
// 5. VISUALISATION — points over the strata-defining map, for a sanity look
// =============================================================================
var roi = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');
Map.centerObject(roi, 9);
Map.setOptions('SATELLITE');
Map.addLayer(pts, {color: 'cyan'}, 'validation points', false);
