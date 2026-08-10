// =============================================================================
// =============================================================================
// 13: AREA AUDIT — diagnostic only (no assets created, no sampling, no exports
//     that touch the pipeline)
// =============================================================================
// =============================================================================
//
// WHY THIS EXISTS
// ---------------
// Two independent numbers for the same thing disagree by ~3x:
//
//   * downstream R analysis (paper Fig. 2): ~980-1640 ha of mapped retama PER YEAR
//   * 09v_valGenerator.js section 5:        5517.8 ha for stratum S1, which is the
//                                           retama mapped in 2017 *OR* 2025
//
// A union of two years cannot exceed the sum of the two years, so at least one of
// them is wrong. This script measures the area every defensible way and prints
// them side by side so the disagreement can be localised instead of argued about.
//
// THE HYPOTHESIS IT TESTS
// -----------------------
// 08_RF2patchFilter.js writes a PRESENCE-ONLY asset. In that file the expression
//
//     .updateMask( pred.eq(0).or(connected.gte(MIN_PATCH_PIX)) )
//
// looks like it preserves class 0, but `connected` comes from `pred.selfMask()`
// and is therefore masked wherever pred == 0. EE's `.or()` does not
// short-circuit: it INTERSECTS masks. So the surviving pixels are exactly
// "class 1 in a patch of >= MIN_PATCH_PIX", and every class-0 pixel is nodata.
//
// Now, 09v_valGenerator.js:157-167 measures the stratum areas with
// `scale: AREA_SCALE` = 100 m and no `crs`/`reproject`. At 100 m, EE reads the
// asset from its image pyramid. Because the asset holds ONLY 1s among valid
// pixels, the pyramid average over valid pixels is 1 in every 100 m cell that
// contains at least one 10 m retama pixel -- and `ee.Image.pixelArea()` at that
// scale returns ~10,000 m2. Net effect: ONE 10 m RETAMA PIXEL PAINTS A WHOLE
// HECTARE.
//
// Section 4 is the decisive test: if the hypothesis holds, the S1 area grows
// monotonically as the measuring scale coarsens. If it does NOT grow, the
// hypothesis is wrong and section 3 still gives the authoritative per-year area.
//
// A SECOND, SMALLER EFFECT
// ------------------------
// Every Export.image.toAsset in this pipeline uses `scale: 10` with NO `crs`, so
// the assets live in EPSG:4326 with square *degree* pixels: ~10 m north-south but
// 10*cos(41 deg) ~ 7.5 m east-west, i.e. ~75.5 m2 of ground, not 100 m2. Counting
// pixels and multiplying by 100 m2 -- which is literally what Methods 2.3 of the
// manuscript describes -- overestimates by ~32%. Section 3 reports both
// conventions so it is clear which one any given number follows.
//
// WHAT THIS SCRIPT DOES NOT DO
// ----------------------------
// It creates no assets, draws no samples and modifies nothing. It is safe to run
// at any time. The validation points, their strata and the human labels are NOT
// affected by the bug above (the stratifiedSample runs at scale 10 in UTM19S);
// only the area WEIGHTS Wi are.
// =============================================================================


// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER    = 'Retamap/Retamap_GEE_Exports';

// Interactive prints can time out on the 10 m full-ROI reductions. Flip this to
// true to also run the whole audit as a batch task (no time limit) that writes a
// tidy CSV to Drive. The printed and exported numbers are identical.
var EXPORT_TO_DRIVE = false;

var RUN_SUFFIX = '_trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k';

var VAL_YEARS = [2017, 2025];   // the two years the 09v strata were built from

// Section 3 runs one full-ROI reduction per year per convention, so it defaults
// to the two validated years -- enough to diagnose, and 2017 is directly
// comparable to the ~1370 ha of Fig. 2A. Widen to all nine years for the full
// Fig. 2 comparison, ideally with EXPORT_TO_DRIVE = true so it runs as a batch
// task instead of against the Code Editor's interactive time limit:
//   var YEARS = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var YEARS = VAL_YEARS;

// Scale ladder for section 4. 100 is the value 09v_valGenerator.js uses today;
// 10 is the native resolution of the maps.
var SCALES = [10, 20, 50, 100, 200, 500];

// Section 7 also checks the 1 km aggregation path. The grid-FC variant needs
// 12_grid_1000m_clipped to exist (12_gridExport.js currently defaults to
// CELL_SIZE = 3000); the raster variant needs no asset.
var CHECK_GRID_FC = false;
var GRID_FC_ASSET = ASSET_PREFIX + '12_grid_1000m_clipped';

var UTM19S  = ee.Projection('EPSG:32719');
var roi     = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');
var roiGeom = roi.geometry();

// Values currently hard-coded in 09v_valAnalysis_compute.R:79-86, transcribed by
// hand from the GEE console on 2026-08-07. These are the numbers under audit.
var OLD_AREA_M2 = {s1: 55178037, s2: 591685230, s3: 9273053977};

// Rows accumulated across all sections -> one print, one optional CSV.
var auditRows = [];
var row = function(section, metric, value) {
  auditRows.push(ee.Feature(null, {
    section: section,
    metric : metric,
    value  : value
  }));
  return value;
};


// =============================================================================
// 1. THE MAPS, AND THE GRID THEY LIVE ON
// =============================================================================
// Mirrors 09v_valGenerator.js:77-83 exactly.
var finalMap = function(year) {
  return ee.Image(ASSET_PREFIX + '08_RF2_prediction_' + year + RUN_SUFFIX)
    .select('classification');
};

var map2017 = finalMap(VAL_YEARS[0]);
var map2025 = finalMap(VAL_YEARS[1]);

// projA = "grid A", the native grid the maps were classified and stored on.
var projA = map2017.projection();

print('--- 1. GRIDS -------------------------------------------------------');
print('08_RF2_prediction_2017 projection:', projA);
print('08_RF2_prediction_2017 nominalScale (m):', projA.nominalScale());
print('01_MergedBands_2017 projection:',
      ee.Image(ASSET_PREFIX + '01_MergedBands_' + VAL_YEARS[0]).projection());

// True ground area of one native pixel, sampled at the ROI centroid. If the
// assets really are EPSG:4326 at "scale 10", this comes out near 75 m2, not 100.
var pixAreaHere = ee.Image.pixelArea().reduceRegion({
  reducer  : ee.Reducer.first(),
  geometry : roiGeom.centroid(10),
  crs      : projA,
  maxPixels: 1e6
}).get('area');
print('Ground area of ONE native pixel at the ROI centroid (m2):', pixAreaHere);
row(1, 'native_pixel_area_m2_at_roi_centroid', pixAreaHere);

// Empirical check of the presence-only claim: for a presence-only asset the
// valid-data footprint and the "== 1" footprint are the SAME. If class 0 were
// stored, the mask footprint would cover the whole ROI instead.
var maskedFootprint = ee.Image.pixelArea().updateMask(map2017.mask())
  .reduceRegion({reducer: ee.Reducer.sum(), geometry: roiGeom, crs: projA,
                 maxPixels: 1e13, tileScale: 4}).get('area');
var onesFootprint = ee.Image.pixelArea().updateMask(map2017.unmask(0).eq(1))
  .reduceRegion({reducer: ee.Reducer.sum(), geometry: roiGeom, crs: projA,
                 maxPixels: 1e13, tileScale: 4}).get('area');
print('Presence-only check 2017 -- valid-data area vs "==1" area (m2). ' +
      'Equal => the asset stores ONLY retama; class 0 is nodata:',
      ee.Dictionary({valid_data: maskedFootprint, equals_one: onesFootprint}));
row(1, 'map2017_valid_data_area_m2', maskedFootprint);
row(1, 'map2017_equals_one_area_m2', onesFootprint);


// =============================================================================
// 2. THE TWO DENOMINATORS
// =============================================================================
// The manuscript quotes a 2,469,653 ha study area (the ROI polygon). The
// validation divides by ~991,992 ha (the valid-observation footprint, after the
// water and elevation masks of 01_sentinelMosaic.js). They are not the same
// thing and a percentage is meaningless until you say which one it is over.
print('--- 2. DENOMINATORS ------------------------------------------------');

var roiPolyArea = roiGeom.area(1);
print('ROI polygon area (ha) -- compare with the 2,469,653 ha in the paper:',
      roiPolyArea.divide(1e4));
row(2, 'roi_polygon_area_ha', roiPolyArea.divide(1e4));

// Mirrors 09v_valGenerator.js:98-100.
var landMask = ee.Image(ASSET_PREFIX + '01_MergedBands_' + VAL_YEARS[0]).select(0).mask()
  .or(ee.Image(ASSET_PREFIX + '01_MergedBands_' + VAL_YEARS[1]).select(0).mask());
var dataMask = landMask;

var dataMaskAreaAt = function(s) {
  return ee.Number(ee.Image.pixelArea().updateMask(dataMask).reduceRegion({
    reducer  : ee.Reducer.sum(),
    geometry : roiGeom,
    scale    : s,
    maxPixels: 1e13,
    tileScale: 4
  }).get('area')).divide(1e4);
};
var dmDict = {};
SCALES.forEach(function(s) {
  dmDict['scale_' + s] = row(2, 'dataMask_area_ha_at_scale_' + s, dataMaskAreaAt(s));
});
print('Valid-observation footprint (ha) by measuring scale. This is the ' +
      'denominator the validation uses (TOTAL_HA = 991,992 today):',
      ee.Dictionary(dmDict));


// =============================================================================
// 3. MAPPED RETAMA AREA PER YEAR -- three conventions, native resolution
// =============================================================================
// (a) pixel count x 100 m2  -- what Methods 2.3 of the manuscript describes
// (b) sum(pixelArea) on the native grid -- the true ground area
// (c) sum(pixelArea) reprojected to UTM19S at 10 m -- sanity check on (b)
//
// (b) / (a) should land near cos(latitude) ~ 0.75 if the EPSG:4326 pixel-shape
// effect is the only difference.
print('--- 3. MAPPED AREA PER YEAR (native resolution) --------------------');

var retamaMask = function(year) {
  return finalMap(year).unmask(0).eq(1);
};

// Convention (a): count pixels, multiply by the nominal 100 m2. Counting the
// mask itself (rather than a constant image) keeps the reduction pinned to the
// asset's own projection instead of the constant's default one.
var areaCountNominal = function(mask) {
  return ee.Number(mask.selfMask().rename('n')
    .reduceRegion({
      reducer  : ee.Reducer.count(),
      geometry : roiGeom,
      crs      : projA,
      maxPixels: 1e13,
      tileScale: 4
    }).get('n')).multiply(100).divide(1e4);
};

// Conventions (b) and (c): true pixel area, on a caller-supplied grid.
// Passing `crs` WITHOUT `scale` makes the reduction run on that projection's own
// transform -- i.e. exactly the grid the asset is stored on, with no resampling.
var areaTrueOn = function(mask, proj) {
  return ee.Number(ee.Image.pixelArea().updateMask(mask).reduceRegion({
    reducer  : ee.Reducer.sum(),
    geometry : roiGeom,
    crs      : proj,
    maxPixels: 1e13,
    tileScale: 4
  }).get('area')).divide(1e4);
};

var perYear = {};
YEARS.forEach(function(y) {
  var m  = retamaMask(y);
  var ha = ee.Dictionary({
    a_count_x100m2 : row(3, 'area_ha_' + y + '_count_x100m2', areaCountNominal(m)),
    b_true_gridA   : row(3, 'area_ha_' + y + '_true_gridA',   areaTrueOn(m, projA)),
    c_true_utm19s  : row(3, 'area_ha_' + y + '_true_utm19s',
                          areaTrueOn(m, UTM19S.atScale(10)))
  });
  perYear['y' + y] = ha;
});
print('Mapped retama area (ha) per year, three conventions. Compare with the ' +
      'paper: 2017 ~1370 ha, 2025 ~1080 ha (Fig. 2A):',
      ee.Dictionary(perYear));


// =============================================================================
// 4. THE SCALE LADDER -- the decisive test
// =============================================================================
// isRetama is built EXACTLY as in 09v_valGenerator.js:88-91, and the area is
// measured with the generator's own code path, varying only `scale`. If the
// pyramid hypothesis is right, S1 inflates monotonically with coarser scale and
// the value at 100 m is several times the value at 10 m.
print('--- 4. SCALE LADDER (the test) -------------------------------------');

var m17b = map2017.unmask(0);
var m25b = map2025.unmask(0);
var isRetama = m17b.eq(1).or(m25b.eq(1));

// Mirrors 09v_valGenerator.js:114-116.
var isNear = isRetama
  .focal_max({radius: 500, units: 'meters', kernelType: 'square'})
  .rename('near');

// Mirrors 09v_valGenerator.js:131-137.
var buildStratum = function(retamaImg, nearImg) {
  var r  = retamaImg.rename('v');
  var n  = nearImg.rename('v');
  var nr = r.multiply(-1).add(1);                          // 1 - r
  return r.add(nr.multiply(n.multiply(-1).add(3)))         // r + (1-r)*(3-n)
    .toByte()
    .updateMask(dataMask)
    .rename('stratum');
};
var stratum = buildStratum(isRetama, isNear);

// The generator's areaOfStratum, verbatim except that `scale` is a parameter.
var areaOfStratumAtScale = function(s, sc) {
  return ee.Number(ee.Image.pixelArea()
    .updateMask(stratum.eq(s))
    .reduceRegion({
      reducer  : ee.Reducer.sum(),
      geometry : roiGeom,
      scale    : sc,
      maxPixels: 1e13,
      tileScale: 4
    }).get('area'));
};

var ladder = {};
SCALES.forEach(function(sc) {
  ladder['scale_' + sc] = ee.Dictionary({
    S1_ha: row(4, 'S1_ha_at_scale_' + sc, areaOfStratumAtScale(1, sc).divide(1e4)),
    S2_ha: row(4, 'S2_ha_at_scale_' + sc, areaOfStratumAtScale(2, sc).divide(1e4)),
    S3_ha: row(4, 'S3_ha_at_scale_' + sc, areaOfStratumAtScale(3, sc).divide(1e4))
  });
});
print('Stratum areas (ha) vs measuring scale. Today the generator reports the ' +
      'scale_100 row (S1 = 5518 ha). Monotonic growth of S1 with scale ' +
      'CONFIRMS the pyramid hypothesis:',
      ee.Dictionary(ladder));

// Inflation factor of the current setting relative to native resolution.
var s1At10  = areaOfStratumAtScale(1, 10);
var s1At100 = areaOfStratumAtScale(1, 100);
print('S1 inflation factor, scale 100 / scale 10 (1.0 would refute the ' +
      'hypothesis):', s1At100.divide(s1At10));
row(4, 'S1_inflation_factor_100_over_10', s1At100.divide(s1At10));


// =============================================================================
// 5. RECOMPUTED Wi -- the numbers that go into the R
// =============================================================================
// The points were drawn with `scale: 10, projection: UTM19S`
// (09v_valGenerator.js:304-305). For the Olofsson estimator the weights must
// describe the SAME strata the points were drawn from, so the UTM19S @10 m row
// below is the one to trust. The grid-A row is reported to show the two agree up
// to the pixel-shape effect.
print('--- 5. RECOMPUTED Wi -----------------------------------------------');

var weightsOn = function(proj, tag) {
  var a1 = areaTrueOn(stratum.eq(1), proj).multiply(1e4);   // back to m2
  var a2 = areaTrueOn(stratum.eq(2), proj).multiply(1e4);
  var a3 = areaTrueOn(stratum.eq(3), proj).multiply(1e4);
  var aT = a1.add(a2).add(a3);
  row(5, tag + '_S1_area_m2', a1);
  row(5, tag + '_S2_area_m2', a2);
  row(5, tag + '_S3_area_m2', a3);
  row(5, tag + '_W1', a1.divide(aT));
  row(5, tag + '_W2', a2.divide(aT));
  row(5, tag + '_W3', a3.divide(aT));
  return ee.Dictionary({
    S1_m2: a1, S2_m2: a2, S3_m2: a3,
    S1_ha: a1.divide(1e4), S2_ha: a2.divide(1e4), S3_ha: a3.divide(1e4),
    W1: a1.divide(aT), W2: a2.divide(aT), W3: a3.divide(aT),
    sumW: a1.add(a2).add(a3).divide(aT),
    // ratio to the values currently baked into 09v_valAnalysis_compute.R
    S1_ratio_vs_current: a1.divide(OLD_AREA_M2.s1),
    S2_ratio_vs_current: a2.divide(OLD_AREA_M2.s2),
    S3_ratio_vs_current: a3.divide(OLD_AREA_M2.s3)
  });
};

print('>> PASTE THESE into 09v_valAnalysis_compute.R section 0 ' +
      '(UTM19S @10 m -- the grid the points were drawn on):',
      weightsOn(UTM19S.atScale(10), 'utm19s10'));
print('Same on the native grid A (cross-check; should differ only by the ' +
      'EPSG:4326 pixel-shape factor):',
      weightsOn(projA, 'gridA'));
print('Currently in the R (measured at scale 100) -- m2:', OLD_AREA_M2);


// =============================================================================
// 6. INTERNAL CONSISTENCY CHECKS
// =============================================================================
// These do not depend on the hypothesis being right. If they fail, the problem
// is in how `stratum` is assembled, not in the measuring scale.
print('--- 6. CONSISTENCY CHECKS ------------------------------------------');

var proj10 = UTM19S.atScale(10);

var a2017  = areaTrueOn(retamaMask(VAL_YEARS[0]), proj10);
var a2025  = areaTrueOn(retamaMask(VAL_YEARS[1]), proj10);
var aUnion = areaTrueOn(retamaMask(VAL_YEARS[0]).or(retamaMask(VAL_YEARS[1])), proj10);
var aInter = areaTrueOn(retamaMask(VAL_YEARS[0]).and(retamaMask(VAL_YEARS[1])), proj10);
var aS1    = areaTrueOn(stratum.eq(1), proj10);

print('Check A -- S1 must equal area(2017 union 2025). Difference in ha ' +
      '(0 expected; non-zero means the dataMask clips part of S1):',
      aS1.subtract(aUnion));
print('Check B -- inclusion-exclusion, A(2017)+A(2025)-A(inter)-A(union). ' +
      'Difference in ha (0 expected):',
      a2017.add(a2025).subtract(aInter).subtract(aUnion));
print('Areas (ha) 2017 / 2025 / union / intersection. The intersection is the ' +
      'part of the map stable across the series endpoints:',
      ee.Dictionary({y2017: a2017, y2025: a2025, union: aUnion, inter: aInter}));

row(6, 'check_S1_minus_union_ha', aS1.subtract(aUnion));
row(6, 'check_inclusion_exclusion_ha',
       a2017.add(a2025).subtract(aInter).subtract(aUnion));
row(6, 'area_ha_union_2017_2025', aUnion);
row(6, 'area_ha_inter_2017_2025', aInter);


// =============================================================================
// 7. THE 1 km AGGREGATION PATH
// =============================================================================
// The downstream R analysis works on a 1000 x 1000 m grid. Aggregating a sparse
// 10 m class to 1 km is exactly where area silently appears or disappears, so
// both the correct and the naive route are measured here. If the naive route
// reproduces the paper's ~1370 ha for 2017, the problem is downstream; if the
// native number in section 3 already is ~1370 ha, then the downstream analysis
// is right and the only thing broken is the generator.
print('--- 7. 1 km AGGREGATION PATH ---------------------------------------');

var oneKm = function(year) {
  var m = retamaMask(year);

  // Correct route: average the fraction of each 1 km cell that is retama, then
  // weight by the cell's area. Area-preserving by construction.
  // maxPixels is generous on purpose: 10 m -> 1 km is 10,000 input pixels
  // exactly, and grid misalignment between EPSG:4326 and UTM pushes it over.
  // Too low a value makes EE fall back to a coarser pyramid level, which is the
  // very failure mode this section is meant to expose.
  var frac = m.reduceResolution({reducer: ee.Reducer.mean(), maxPixels: 65535})
    .reproject({crs: UTM19S, scale: 1000});
  var haCorrect = ee.Number(frac.multiply(ee.Image.pixelArea()).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: roiGeom, crs: UTM19S, scale: 1000,
    maxPixels: 1e13, tileScale: 4
  }).get('area')).divide(1e4);

  // Naive route: reproject straight to 1 km with no reduceResolution, so EE
  // resamples from the pyramid. This is the classic way to lose (or invent) a
  // sparse class wholesale.
  var haNaive = ee.Number(ee.Image.pixelArea()
    .updateMask(m.reproject({crs: UTM19S, scale: 1000}))
    .reduceRegion({
      reducer: ee.Reducer.sum(), geometry: roiGeom, crs: UTM19S, scale: 1000,
      maxPixels: 1e13, tileScale: 4
    }).get('area')).divide(1e4);

  row(7, 'area_ha_' + year + '_1km_reduceResolution', haCorrect);
  row(7, 'area_ha_' + year + '_1km_naive_reproject', haNaive);
  return ee.Dictionary({reduceResolution: haCorrect, naive_reproject: haNaive});
};

var kmDict = {};
VAL_YEARS.forEach(function(y) { kmDict['y' + y] = oneKm(y); });
print('Mapped area (ha) via the 1 km grid, correct vs naive aggregation:',
      ee.Dictionary(kmDict));

if (CHECK_GRID_FC) {
  // Same thing through the actual grid FeatureCollection the R analysis uses.
  var grid = ee.FeatureCollection(GRID_FC_ASSET);
  VAL_YEARS.forEach(function(y) {
    var perCell = ee.Image.pixelArea().updateMask(retamaMask(y))
      .reduceRegions({
        collection: grid,
        reducer   : ee.Reducer.sum(),
        crs       : proj10,
        tileScale : 4
      });
    var tot = ee.Number(perCell.aggregate_sum('sum')).divide(1e4);
    print('Mapped area (ha) ' + y + ' summed over ' + GRID_FC_ASSET + ':', tot);
    row(7, 'area_ha_' + y + '_gridFC_sum', tot);
  });
}


// =============================================================================
// 8. TIDY OUTPUT
// =============================================================================
var auditFc = ee.FeatureCollection(auditRows);
print('--- 8. ALL ROWS ----------------------------------------------------');
print('Full audit table (section, metric, value):', auditFc);

if (EXPORT_TO_DRIVE) {
  Export.table.toDrive({
    collection    : auditFc,
    description   : 'Export_13_areaAudit',
    folder        : DRIVE_FOLDER,
    fileNamePrefix: '13_areaAudit' + RUN_SUFFIX,
    fileFormat    : 'CSV',
    selectors     : ['section', 'metric', 'value']
  });
}
