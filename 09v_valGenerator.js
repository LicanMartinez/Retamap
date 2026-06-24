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

var UTM19S  = ee.Projection('EPSG:32719');
var roi     = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');
var roiGeom = roi.geometry();

// ── Olofsson sample-size parameters (editable) ───────────────────────────────
// n_i = U_i (1 - U_i) / SE^2 , with a per-stratum floor. U_i = expected user's
// accuracy of the mapped class. These are the SUGGESTED sizes; you may override
// N_OVERRIDE below (e.g. raise S2 to hunt false positives with more power).
var EXP_UA    = {s1: 0.85, s2: 0.85, s3: 0.9};
var TARGET_SE = 0.02;
var FLOOR     = {s1: 0, s2: 0, s3: 0};
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
  s1: nFormula(EXP_UA.s1, TARGET_SE, FLOOR.s1),
  s2: nFormula(EXP_UA.s2, TARGET_SE, FLOOR.s2),
  s3: nFormula(EXP_UA.s3, TARGET_SE, FLOOR.s3)
};
var nFinal = {
  s1: (N_OVERRIDE.s1 === null) ? nSuggest.s1 : N_OVERRIDE.s1,
  s2: (N_OVERRIDE.s2 === null) ? nSuggest.s2 : N_OVERRIDE.s2,
  s3: (N_OVERRIDE.s3 === null) ? nSuggest.s3 : N_OVERRIDE.s3
};
print('Olofsson per-stratum n (formula):', nSuggest);
print('Per-stratum n USED (after override):', nFinal,
      'TOTAL =', nFinal.s1 + nFinal.s2 + nFinal.s3);

// Olofsson Eq.13 cross-check: total n to reach a target SE of OVERALL accuracy.
// n_total = ( Σ Wi·Si / S(OA) )^2 , Si = sqrt(Ui(1-Ui)).
var Si = function(u) { return ee.Number(u * (1 - u)).sqrt(); };
var sumWS = ee.Number(W.get('s1')).multiply(Si(EXP_UA.s1))
  .add(ee.Number(W.get('s2')).multiply(Si(EXP_UA.s2)))
  .add(ee.Number(W.get('s3')).multiply(Si(EXP_UA.s3)));
print('Olofsson n_total for S(OA)=' + TARGET_SE + ' (cross-check):',
      sumWS.divide(TARGET_SE).pow(2).ceil());

// =============================================================================
// 7. STRATIFIED SAMPLE (pixel centroids) + carry labels + lon/lat
// =============================================================================
var samples = stratStack.stratifiedSample({
  numPoints  : 0,                       // ignored when classPoints is given
  classBand  : 'stratum',
  region     : roiGeom,
  scale      : 10,
  projection : UTM19S,
  classValues: [1, 2, 3],
  classPoints: [nFinal.s1, nFinal.s2, nFinal.s3],
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
// 9. EXPORT raw points (WITH labels) → Drive CSV for 09v_valRandomize.R
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
