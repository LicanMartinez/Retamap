// =============================================================================
// =============================================================================
// 09vNear: Independent-validation point GENERATOR — NEAR-ONLY variant (2 strata)
// =============================================================================
// =============================================================================
//
// Variant of 09v_valGenerator.js. Same reference-label strategy (independent
// human photo-interpretation, sampling stratified BY THE FINAL MAP Stage C,
// Olofsson et al. 2014), but DROPS stratum S3 (far non-retama background).
//
// Rationale: outside the ~1 km neighbourhood of mapped retama (S2), we trust
// the model's commission error is ~0 — false positives only show up near
// roads/urbanization/agriculture adjacent to actual retama patches, not in the
// bulk far background. Validating S3 spends most of the point budget on a
// stratum whose error rate is already known to be negligible, and the 3-way
// Wi area-weighting (S3 is by far the largest stratum) makes UA/PA harder to
// read as plain commission/omission rates. Dropping it also makes point
// collection far easier (no need to scatter validators across the whole ROI).
//
// CONSEQUENCE — read this before interpreting results: dropping S3 shrinks the
// assessed POPULATION from "the whole ROI" to "the near-retama domain"
// (S1 ∪ S2, i.e. mapped retama + its ~1 km buffer). All Wi, OA, UA, PA, and
// area estimates produced by this script/its R companion are CONDITIONAL on
// that restricted domain — they are NOT an unbiased estimator for the full
// study area. In particular, PA (omission) only "sees" omission errors that
// happen near already-mapped retama; a hypothetical retama patch missed
// entirely (isolated, no neighbour pixel mapped as retama) would fall in S3
// and is invisible to this variant. That blind spot is the explicit trade-off
// being made here (see rationale above) — use 09v_valGenerator.js (3-strata,
// full-ROI) if/when that blind spot needs to be checked.
//
// Two strata:
//   S1 — mapped retama : pixel = 1 in 2017 AND/OR 2025
//   S2 — near non-retama: not-S1, within ~1 km of S1 (focal_max dilation)
//
// Output: Drive CSV '09vNear_valPoints_raw{RUN_SUFFIX}.csv', same row shape as
// the 3-strata generator (stratum, m2017, m2025, lon, lat), consumed by
// 09vNear_valRandomize.R.
// =============================================================================


// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER    = 'Retamap/Retamap_GEE_Exports';
var EXPORT_TO_DRIVE = true;   // generator's whole purpose is the Drive CSV

var RUN_SUFFIX = '_trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k';

var VAL_YEARS = [2017, 2025];   // years validated; also the years that define S2

var SEED       = 42;
var AREA_SCALE = 100;   // m, coarse scale for stratum-area (Wi) estimation only

var UTM19S  = ee.Projection('EPSG:32719');
var roi     = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');
var roiGeom = roi.geometry();

// ── Olofsson sample-size parameters (editable) ───────────────────────────────
// n_i = U_i (1 - U_i) / SE_i^2 , per-stratum floor. Only s1/s2 now — same
// defaults as the 3-strata generator's S1/S2 entries.
var EXP_UA    = {s1: 0.85, s2: 0.85};
var TARGET_SE = {s1: 0.02, s2: 0.04};
var FLOOR     = {s1: 0, s2: 0};
// Set any entry to a number to override the formula; null → use formula value.
var N_OVERRIDE = {s1: null, s2: null};

// =============================================================================
// 1. LOAD FINAL MAPS (Stage C)
// =============================================================================
var finalMap = function(year) {
  return ee.Image(ASSET_PREFIX + '08_RF2_prediction_' + year + RUN_SUFFIX)
    .select('classification');
};

var map2017 = finalMap(2017);
var map2025 = finalMap(2025);

// =============================================================================
// 2. S1 CONDITION — retama in 2017 and/or 2025
// =============================================================================
var m17b = map2017.unmask(0);
var m25b = map2025.unmask(0);

var isRetama = m17b.eq(1).or(m25b.eq(1));            // S1 condition; reused below

// Data mask: the 08_RF2_prediction assets are PRESENCE maps — only retama
// (value 1) is stored, everything else is nodata (NOT value 0). So the asset
// cannot define the "valid territory" for the non-retama stratum. Instead,
// take the valid-land mask (water + elevation already excluded) from
// 01_MergedBands, which is the actual observation footprint the model ran on.
var landMask = ee.Image(ASSET_PREFIX + '01_MergedBands_2017').select(0).mask()
  .or(ee.Image(ASSET_PREFIX + '01_MergedBands_2025').select(0).mask());
var dataMask = landMask;

// =============================================================================
// 3. NEAR-RETAMA MASK — within ~1 km of retama; ALSO defines the domain
// =============================================================================
// focal_max dilates isRetama by 500 m in each direction (≈ 1 km square window).
// isNear already contains isRetama itself (dilation preserves the seed pixels),
// so isNear alone is the full S1∪S2 domain — no separate "far" stratum exists
// in this variant; pixels outside isNear are simply masked out (not sampled).
var isNear = isRetama
  .focal_max({radius: 500, units: 'meters', kernelType: 'square'})
  .rename('near');

// =============================================================================
// 4. BUILD STRATUM IMAGE {1,2} + carry both map labels
// =============================================================================
// Arithmetic only (avoids .where()/.and() projection-ambiguity issues seen in
// the 3-strata script): retama → 1; near non-retama → 2.
//   stratum = 2 − isRetama
var stratum = ee.Image(2).subtract(isRetama)
  .toByte()
  .updateMask(isNear)     // restrict domain to S1 ∪ S2 (drops S3 entirely)
  .updateMask(dataMask)
  .rename('stratum');

// Diagnostic: should show {1: <count>, 2: <count>}
// print('Stratum pixel histogram (sanity check):',
//   stratum.reduceRegion({
//     reducer  : ee.Reducer.frequencyHistogram(),
//     geometry : roiGeom,
//     scale    : AREA_SCALE,
//     maxPixels: 1e13,
//     tileScale: 4
//   }));

// Bands carrying the per-year final-map label (masked → property absent → NA).
var stratStack = stratum
  .addBands(map2017.rename('m2017'))
  .addBands(map2025.rename('m2025'));

// =============================================================================
// 5. STRATUM AREAS → WEIGHTS Wi  (coarse scale; weights only)
// =============================================================================
// NOTE: Wi here sum to 1 over the RESTRICTED domain (S1∪S2), not the whole
// ROI — see the CONSEQUENCE note in the header. They are area-within-domain
// weights, used only to combine S1/S2 into OA/PA estimates for that domain.
var areaOfStratum = function(s) {
  return ee.Number(ee.Image.pixelArea()
    .updateMask(stratum.eq(s))
    .reduceRegion({
      reducer  : ee.Reducer.sum(),
      geometry : roiGeom,
      scale    : AREA_SCALE,
      maxPixels: 1e13,
      tileScale: 4
    }).get('area'));
};
var a1 = areaOfStratum(1), a2 = areaOfStratum(2);
print('Stratum areas within the near-domain (m^2):', ee.Dictionary({s1: a1, s2: a2}));
var aTot = a1.add(a2);
var W = ee.Dictionary({s1: a1.divide(aTot), s2: a2.divide(aTot)});
print('Area weights Wi within near-domain (should sum to 1):', W,
      'Σ=', ee.Number(W.get('s1')).add(W.get('s2')));

// =============================================================================
// 6. OLOFSSON SAMPLE SIZES
// =============================================================================
var nFormula = function(u, se, floor) {
  return Math.max(Math.ceil((u * (1 - u)) / (se * se)), floor);
};
var nSuggest = {
  s1: nFormula(EXP_UA.s1, TARGET_SE.s1, FLOOR.s1),
  s2: nFormula(EXP_UA.s2, TARGET_SE.s2, FLOOR.s2)
};
var nFinal = {
  s1: (N_OVERRIDE.s1 === null) ? nSuggest.s1 : N_OVERRIDE.s1,
  s2: (N_OVERRIDE.s2 === null) ? nSuggest.s2 : N_OVERRIDE.s2
};
print('Olofsson per-stratum n (formula):', nSuggest);
print('Per-stratum n USED (after override):', nFinal,
      'TOTAL =', nFinal.s1 + nFinal.s2);

// Cross-check: achieved SE of OVERALL accuracy (within the near-domain) given
// nFinal: SE(OA) = sqrt( Σ Wi² · Ui(1-Ui) / n_i ), 2-term version.
var varTerm = function(s) {
  return ee.Number(W.get(s)).pow(2)
    .multiply(EXP_UA[s] * (1 - EXP_UA[s]))
    .divide(nFinal[s]);
};
var achievedSE_OA = varTerm('s1').add(varTerm('s2')).sqrt();
print('Achieved overall-accuracy SE given nFinal (cross-check, near-domain only):',
      achievedSE_OA);

// =============================================================================
// 6b. EXPECTED SE OF PRODUCER'S ACCURACY (delta method) — 2-stratum version
// =============================================================================
// Same logic as the 3-strata script's §6b, with the S3 term dropped. PA here
// is PA *within the near-domain* (S1∪S2) — see header CONSEQUENCE note.
//
// q_i = P(truly retama | stratum i):  S1: q = UA_S1   S2: q = 1−UA_S2 (omission)
// PA_ret = W1·q1 / (W1·q1 + W2·q2)
// V(PA)  = Σ (∂PA/∂q_i)² · q_i(1−q_i)/n_i

var _W1 = ee.Number(W.get('s1'));
var _W2 = ee.Number(W.get('s2'));

// ── PA retama ────────────────────────────────────────────────────────────────
var qR1 = EXP_UA.s1,  qR2 = 1 - EXP_UA.s2;

var paN  = _W1.multiply(qR1);
var paD  = paN.add(_W2.multiply(qR2));
var paRet = paN.divide(paD);

var paD2 = paD.pow(2);
var dr1 = _W1.multiply(paD.subtract(paN)).divide(paD2);     // ∂PA/∂q1
var dr2 = paN.multiply(-1).multiply(_W2).divide(paD2);      // ∂PA/∂q2

var vr1 = qR1*(1-qR1)/nFinal.s1;
var vr2 = qR2*(1-qR2)/nFinal.s2;

var cR1 = dr1.pow(2).multiply(vr1);
var cR2 = dr2.pow(2).multiply(vr2);
var vPaRet = cR1.add(cR2);
var sePaRet = vPaRet.sqrt();

print('─── Expected PA & SE (delta method, near-domain only) ───');
print('RETAMA  — E[PA]:', paRet, '  SE(PA):', sePaRet,
      '  (cf. SE(UA) = ' + TARGET_SE.s1 + ')');
print('  V(PA) % by stratum →  S1:', cR1.divide(vPaRet).multiply(100),
      '  S2:', cR2.divide(vPaRet).multiply(100));

// ── PA non-retama ────────────────────────────────────────────────────────────
var qN1 = 1 - EXP_UA.s1,  qN2 = EXP_UA.s2;
var paNn = _W2.multiply(qN2);
var paNd = _W1.multiply(qN1).add(paNn);
var paNR = paNn.divide(paNd);

var paNd2 = paNd.pow(2);
var dn1 = paNn.multiply(-1).multiply(_W1).divide(paNd2);    // ∂PA/∂q10
var dn2 = _W2.multiply(_W1.multiply(qN1)).divide(paNd2);    // ∂PA/∂q20

var vn1 = qN1*(1-qN1)/nFinal.s1;
var vn2 = qN2*(1-qN2)/nFinal.s2;

var cN1 = dn1.pow(2).multiply(vn1);
var cN2 = dn2.pow(2).multiply(vn2);
var vPaNR = cN1.add(cN2);
var sePaNR = vPaNR.sqrt();

print('NON-RET — E[PA]:', paNR, '  SE(PA):', sePaNR);
print('  V(PA) % by stratum →  S1:', cN1.divide(vPaNR).multiply(100),
      '  S2:', cN2.divide(vPaNR).multiply(100));

// =============================================================================
// 7. STRATIFIED SAMPLE (pixel centroids) + carry labels + lon/lat
// =============================================================================
var samples = stratStack.stratifiedSample({
  numPoints  : 0,                       // ignored when classPoints is given
  classBand  : 'stratum',
  region     : roiGeom,
  scale      : 10,
  projection : UTM19S,
  classValues: [1, 2],
  classPoints: [nFinal.s1, nFinal.s2],
  seed       : SEED,
  geometries : true,
  dropNulls  : false,   // keep points whose label is NA in one year (masked)
  tileScale  : 4
});

// Append lon/lat (centroid coordinates) as plain columns for the R step.
samples = samples.map(function(f) {
  var c = ee.List(f.geometry().transform('EPSG:4326', 1).coordinates());
  return f.set('lon', c.get(0), 'lat', c.get(1));
});

// print('Sampled points per stratum:', samples.aggregate_histogram('stratum'));
// print('First sampled features:', samples.limit(5));

// =============================================================================
// 8. VISUALISATION
// =============================================================================
Map.centerObject(roi, 9);
Map.setOptions('SATELLITE');
Map.addLayer(stratum, {min: 1, max: 2, palette: ['e31a1c', 'ff7f00']},
             'Strata (1=retama,2=near)', true, 0.6);
var palByStratum = function(s, color, name) {
  Map.addLayer(samples.filter(ee.Filter.eq('stratum', s)),
               {color: color}, name, true);
};
palByStratum(1, 'red',    'pts S1 retama');
palByStratum(2, 'orange', 'pts S2 near');

// =============================================================================
// 9. EXPORT raw points (WITH labels) → Drive CSV for 09vNear_valRandomize.R
// =============================================================================
if (EXPORT_TO_DRIVE) {
  Export.table.toDrive({
    collection    : samples,
    description   : '09vNear_valPoints_raw' + RUN_SUFFIX,
    folder        : DRIVE_FOLDER,
    fileNamePrefix: '09vNear_valPoints_raw' + RUN_SUFFIX,
    fileFormat    : 'CSV',
    selectors     : ['stratum', 'm2017', 'm2025', 'lon', 'lat']
  });
}
