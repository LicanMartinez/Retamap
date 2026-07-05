// =============================================================================
// =============================================================================
// 09v: Independent-validation point GENERATOR (map-stratified sampling)
// =============================================================================
// =============================================================================
//
// New validation strategy (replaces the internal holdout of 09_holdoutValidation,
// now archived in experimental/). Reference labels come from independent human
// photo-interpretation; the sampling design is stratified BY THE FINAL MAP
// (Stage C, 08_RF2_prediction) following Olofsson et al. 2014.
//
// Validated years: 2017 and 2025 (start / end of the series). A single point set
// is drawn for both years jointly; each validator labels each year separately.
//
// Three strata that PARTITION the territory (enables an unbiased Olofsson area
// estimator at the scale of the whole study area):
//   S1 — mapped retama : pixel = 1 in 2017 AND/OR 2025
//   S2 — near non-retama: not-S1, INSIDE 1 km cells (via reduceResolution of
//                         isRetama to 1 km) that contain retama in EITHER
//                         validated year (2017 and/or 2025) → FP-prone zones
//   S3 — far non-retama : not-S1, OUTSIDE those cells (bulk of the area)
//
// Output: Drive CSV '09v_valPoints_raw{RUN_SUFFIX}.csv' with one row per sampled
// pixel CENTROID, carrying stratum + both map labels (m2017, m2025) + lon/lat.
// This raw CSV (WITH labels) is consumed locally by 09v_valRandomize.R, which
// strips the labels before building the blind asset FeatureCollection.
// =============================================================================


// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER    = 'Retamap/Retamap_GEE_Exports';
var EXPORT_TO_DRIVE = true;   // generator's whole purpose is the Drive CSV

var RUN_SUFFIX = '_trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k';

var VAL_YEARS = [2017, 2025];   // years validated; also the years that define S2

// No grid-FC prerequisite needed: the near-retama mask is built by
// reduceResolution (raster-only), not by rasterizing the grid FC.

var SEED       = 42;
var AREA_SCALE = 100;   // m, coarse scale for stratum-area (Wi) estimation only

// ── Minimum-distance grid (1 point per cell per stratum) ────────────────────
// Quantize the map into GRID_CELL_SIZE cells; after sampling keep at most one
// point per cell PER STRATUM (via distinct('cellId')) → ~GRID_CELL_SIZE spacing
// within each stratum. Approximate: two points in adjacent cells can fall closer
// than GRID_CELL_SIZE near a shared cell edge. Dedup is per-stratum, so an S1 and
// an S2 point may still be <GRID_CELL_SIZE apart (accepted).
var GRID_CELL_SIZE    = 500;   // metres — target spacing within each stratum
var OVERSAMPLE_FACTOR = 5;     // draw N×nFinal, then keep 1 per cell via distinct()

var UTM19S  = ee.Projection('EPSG:32719');
var roi     = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');
var roiGeom = roi.geometry();

// ── Olofsson sample-size parameters (editable) ───────────────────────────────
// n_i = max( U_i (1 - U_i) / SE_i^2 , FLOOR_i ). U_i = expected user's accuracy of
// the mapped class; SE_i = target SE of that stratum's UA; FLOOR_i = minimum n per
// stratum (safety net). These are the SUGGESTED sizes; you may override via
// N_OVERRIDE below (e.g. raise S2 to hunt false positives with more power). Current
// config: SE = 0.02 on all three strata (tight, equal) + FLOOR = 100 each → formula
// gives ~400/400/119 (total ~919); the floor only binds if a stratum has few 500 m
// cells (plausible for S1 / retama).
var EXP_UA    = {s1: 0.8, s2: 0.8, s3: 0.95};
var TARGET_SE = {s1: 0.02, s2: 0.02, s3: 0.02};
var FLOOR     = {s1: 100, s2: 100, s3: 100};
// Set any entry to a number to override the formula; null → use formula value.
var N_OVERRIDE = {s1: null, s2: null, s3: null};

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
// cannot define the "valid territory" for the non-retama strata. Instead, take
// the valid-land mask (water + elevation already excluded) from 01_MergedBands,
// which is the actual observation footprint the model ran on.
var landMask = ee.Image(ASSET_PREFIX + '01_MergedBands_2017').select(0).mask()
  .or(ee.Image(ASSET_PREFIX + '01_MergedBands_2025').select(0).mask());
var dataMask = landMask;

// print('DEBUG dataMask area (m^2, should be full ROI, >> 55 M):',
//   ee.Image.pixelArea().updateMask(dataMask).reduceRegion({
//     reducer: ee.Reducer.sum(), geometry: roiGeom, scale: AREA_SCALE,
//     maxPixels: 1e13, tileScale: 4
//   }).get('area'));

// =============================================================================
// 3. NEAR-RETAMA MASK — non-retama within ~1 km of retama (FP-prone zone)
// =============================================================================
// focal_max dilates isRetama by 500 m in each direction (≈ 1 km square window).
// Stays entirely in the native Sentinel-2 projection — no vector rasterization,
// no reduceResolution, no reproject.
var isNear = isRetama
  .focal_max({radius: 500, units: 'meters', kernelType: 'square'})
  .rename('near');

// print('isNear area (m^2, sanity check — should be hundreds of km^2):',
//   ee.Image.pixelArea().updateMask(isNear).reduceRegion({
//     reducer: ee.Reducer.sum(), geometry: roiGeom, scale: AREA_SCALE,
//     maxPixels: 1e13, tileScale: 4
//   }).get('area'));

// =============================================================================
// 4. BUILD STRATUM IMAGE {1,2,3} + carry both map labels
// =============================================================================
// Arithmetic combination avoiding .where()/.and()/.not() that had projection
// ambiguity issues. All operands share band name 'v'.
//   retama → 1;  not-retama, near → 2;  not-retama, far → 3
// Formula: stratum = r + (1 − r) × (3 − n)
var r  = isRetama.rename('v');
var n  = isNear.rename('v');
var nr = r.multiply(-1).add(1);              // 1 − r  (1 where not retama)
var stratum = r.add(nr.multiply(n.multiply(-1).add(3)))  // r + (1-r)*(3-n)
  .toByte()
  .updateMask(dataMask)
  .rename('stratum');

// Diagnostic: should show {1: <count>, 2: <count>, 3: <count>}
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
var a1 = areaOfStratum(1), a2 = areaOfStratum(2), a3 = areaOfStratum(3);
print('Stratum areas (m^2):', ee.Dictionary({s1: a1, s2: a2, s3: a3}));
var aTot = a1.add(a2).add(a3);
var W = ee.Dictionary({s1: a1.divide(aTot), s2: a2.divide(aTot), s3: a3.divide(aTot)});
print('Area weights Wi (should sum to 1):', W, 'Σ=', ee.Number(W.get('s1')).add(W.get('s2')).add(W.get('s3')));

// =============================================================================
// 6. OLOFSSON SAMPLE SIZES
// =============================================================================
// Per-stratum n from target SE of that class's user's accuracy (client-side).
var nFormula = function(u, se, floor) {
  return Math.max(Math.ceil((u * (1 - u)) / (se * se)), floor);
};
var nSuggest = {
  s1: nFormula(EXP_UA.s1, TARGET_SE.s1, FLOOR.s1),
  s2: nFormula(EXP_UA.s2, TARGET_SE.s2, FLOOR.s2),
  s3: nFormula(EXP_UA.s3, TARGET_SE.s3, FLOOR.s3)
};
var nFinal = {
  s1: (N_OVERRIDE.s1 === null) ? nSuggest.s1 : N_OVERRIDE.s1,
  s2: (N_OVERRIDE.s2 === null) ? nSuggest.s2 : N_OVERRIDE.s2,
  s3: (N_OVERRIDE.s3 === null) ? nSuggest.s3 : N_OVERRIDE.s3
};
print('Olofsson per-stratum n (formula):', nSuggest);
print('Per-stratum n USED (after override):', nFinal,
      'TOTAL =', nFinal.s1 + nFinal.s2 + nFinal.s3);

// Cross-check: achieved SE of OVERALL accuracy given nFinal (stratified
// variance, Olofsson et al. 2014 Eq.8): SE(OA) = sqrt( Σ Wi² · Ui(1-Ui) / n_i ).
// With 3 independent per-stratum SE targets there's no single target left to
// solve Eq.13 for n_total, so this instead reports the overall-accuracy
// precision that nFinal actually buys (run the formula backward).
var varTerm = function(s) {
  return ee.Number(W.get(s)).pow(2)
    .multiply(EXP_UA[s] * (1 - EXP_UA[s]))
    .divide(nFinal[s]);
};
var achievedSE_OA = varTerm('s1').add(varTerm('s2')).add(varTerm('s3')).sqrt();
print('Achieved overall-accuracy SE given nFinal (cross-check):', achievedSE_OA);

// =============================================================================
// 6b. EXPECTED SE OF PRODUCER'S ACCURACY (delta method)
// =============================================================================
// §6 targets SE of User's Accuracy per stratum. But PA depends on how omission
// errors in OTHER strata combine, weighted by Wi². Delta method on the ratio
// PA_j = p_jj / p_·j  gives SE(PA) for each reference class, so you can tune
// TARGET_SE until PA precision matches UA precision.
//
// Binary reference classes:  retama (1) / non-retama (0)
// q_i = P(truly retama | stratum i):
//   S1: q = UA_S1          S2: q = 1−UA_S2 (omission)   S3: q = 1−UA_S3
//
// PA_ret = W1·q1 / (W1·q1 + W2·q2 + W3·q3)
// V(PA)  = Σ (∂PA/∂q_i)² · q_i(1−q_i)/n_i

var _W1 = ee.Number(W.get('s1'));
var _W2 = ee.Number(W.get('s2'));
var _W3 = ee.Number(W.get('s3'));

// ── PA retama ────────────────────────────────────────────────────────────────
var qR1 = EXP_UA.s1,  qR2 = 1 - EXP_UA.s2,  qR3 = 1 - EXP_UA.s3;

var paN  = _W1.multiply(qR1);
var paD  = paN.add(_W2.multiply(qR2)).add(_W3.multiply(qR3));
var paRet = paN.divide(paD);

var paD2 = paD.pow(2);
var dr1 = _W1.multiply(paD.subtract(paN)).divide(paD2);     // ∂PA/∂q1
var dr2 = paN.multiply(-1).multiply(_W2).divide(paD2);      // ∂PA/∂q2
var dr3 = paN.multiply(-1).multiply(_W3).divide(paD2);      // ∂PA/∂q3

var vr1 = qR1*(1-qR1)/nFinal.s1;
var vr2 = qR2*(1-qR2)/nFinal.s2;
var vr3 = qR3*(1-qR3)/nFinal.s3;

var cR1 = dr1.pow(2).multiply(vr1);
var cR2 = dr2.pow(2).multiply(vr2);
var cR3 = dr3.pow(2).multiply(vr3);
var vPaRet = cR1.add(cR2).add(cR3);
var sePaRet = vPaRet.sqrt();

print('─── Expected PA & SE (delta method) ────────────────────');
print('RETAMA  — E[PA]:', paRet, '  SE(PA):', sePaRet,
      '  (cf. SE(UA) = ' + TARGET_SE.s1 + ')');
print('  V(PA) % by stratum →  S1:', cR1.divide(vPaRet).multiply(100),
      '  S2:', cR2.divide(vPaRet).multiply(100),
      '  S3:', cR3.divide(vPaRet).multiply(100));

// ── PA non-retama ────────────────────────────────────────────────────────────
var qN1 = 1 - EXP_UA.s1,  qN2 = EXP_UA.s2,  qN3 = EXP_UA.s3;
var paNn = _W2.multiply(qN2).add(_W3.multiply(qN3));
var paNd = _W1.multiply(qN1).add(paNn);
var paNR = paNn.divide(paNd);

var paNd2 = paNd.pow(2);
var dn1 = paNn.multiply(-1).multiply(_W1).divide(paNd2);    // ∂PA/∂q10
var dn2 = _W2.multiply(_W1.multiply(qN1)).divide(paNd2);    // ∂PA/∂q20
var dn3 = _W3.multiply(_W1.multiply(qN1)).divide(paNd2);    // ∂PA/∂q30

var vn1 = qN1*(1-qN1)/nFinal.s1;
var vn2 = qN2*(1-qN2)/nFinal.s2;
var vn3 = qN3*(1-qN3)/nFinal.s3;

var cN1 = dn1.pow(2).multiply(vn1);
var cN2 = dn2.pow(2).multiply(vn2);
var cN3 = dn3.pow(2).multiply(vn3);
var vPaNR = cN1.add(cN2).add(cN3);
var sePaNR = vPaNR.sqrt();

print('NON-RET — E[PA]:', paNR, '  SE(PA):', sePaNR);
print('  V(PA) % by stratum →  S1:', cN1.divide(vPaNR).multiply(100),
      '  S2:', cN2.divide(vPaNR).multiply(100),
      '  S3:', cN3.divide(vPaNR).multiply(100));

// =============================================================================
// 7. CELL-ID RASTER (for minimum-distance enforcement)
// =============================================================================
// Quantize UTM coordinates into GRID_CELL_SIZE cells. All 10 m pixels within the
// same cell share a single cellId → after sampling, distinct('cellId') per stratum
// keeps at most 1 point per cell (~GRID_CELL_SIZE spacing). Raster-only, no vector
// grid asset needed.
var coords = ee.Image.pixelCoordinates(UTM19S);
var cellX  = coords.select('x').divide(GRID_CELL_SIZE).floor().toInt32();
var cellY  = coords.select('y').divide(GRID_CELL_SIZE).floor().toInt32();
var cellId = cellX.multiply(100000).add(cellY).rename('cellId');

var stratStackWithCell = stratStack.addBands(cellId);

// =============================================================================
// 8. STRATIFIED SAMPLE — oversample + deduplicate by cell + lon/lat
// =============================================================================
// Step 1: draw OVERSAMPLE_FACTOR × nFinal points per stratum (random placement).
var oversampled = stratStackWithCell.stratifiedSample({
  numPoints  : 0,                       // ignored when classPoints is given
  classBand  : 'stratum',
  region     : roiGeom,
  scale      : 10,
  projection : UTM19S,
  classValues: [1, 2, 3],
  classPoints: [nFinal.s1 * OVERSAMPLE_FACTOR,
                nFinal.s2 * OVERSAMPLE_FACTOR,
                nFinal.s3 * OVERSAMPLE_FACTOR],
  seed       : SEED,
  geometries : true,
  dropNulls  : false,   // keep points whose label is NA in one year (masked)
  tileScale  : 4
});

// Step 2: per stratum, randomize → keep 1 point per cell → take nFinal.
var sampleStratum = function(s, n) {
  return oversampled
    .filter(ee.Filter.eq('stratum', s))
    .randomColumn('rand', SEED + s)
    .sort('rand')
    .distinct('cellId')
    .limit(n);
};

var samples = sampleStratum(1, nFinal.s1)
  .merge(sampleStratum(2, nFinal.s2))
  .merge(sampleStratum(3, nFinal.s3));

// Append lon/lat (centroid coordinates) as plain columns for the R step.
samples = samples.map(function(f) {
  var c = ee.List(f.geometry().transform('EPSG:4326', 1).coordinates());
  return f.set('lon', c.get(0), 'lat', c.get(1));
});

print('Desired n per stratum:', nFinal);
print('Actual n per stratum after min-distance filter:',
      samples.aggregate_histogram('stratum'));

// =============================================================================
// 9. VISUALISATION
// =============================================================================
Map.centerObject(roi, 9);
Map.setOptions('SATELLITE');

// Grid cells (vector, for visual inspection only — logic uses the raster cellId)
var gridVis = roiGeom.bounds(1, UTM19S).coveringGrid(UTM19S, GRID_CELL_SIZE);
Map.addLayer(gridVis, {color: 'aaaaaa', fillColor: '00000000'},
             GRID_CELL_SIZE + ' m grid', false);

Map.addLayer(stratum, {min: 1, max: 3, palette: ['e31a1c', 'ff7f00', '33a02c']},
             'Strata (1=retama,2=near,3=far)', true, 0.6);
Map.addLayer(isNear.selfMask(), {palette: ['1f78b4']},
             'Near-retama zone', false, 0.3);
var palByStratum = function(s, color, name) {
  Map.addLayer(samples.filter(ee.Filter.eq('stratum', s)),
               {color: color}, name, true);
};
palByStratum(1, 'red',    'pts S1 retama');
palByStratum(2, 'orange', 'pts S2 near');
palByStratum(3, 'green',  'pts S3 far');

// =============================================================================
// 10. EXPORT raw points (WITH labels) → Drive CSV for 09v_valRandomize.R
// =============================================================================
if (EXPORT_TO_DRIVE) {
  Export.table.toDrive({
    collection    : samples,
    description   : '09v_valPoints_raw' + RUN_SUFFIX,
    folder        : DRIVE_FOLDER,
    fileNamePrefix: '09v_valPoints_raw' + RUN_SUFFIX,
    fileFormat    : 'CSV',
    selectors     : ['stratum', 'm2017', 'm2025', 'lon', 'lat']
  });
}
