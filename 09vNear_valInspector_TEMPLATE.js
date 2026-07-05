// =============================================================================
// =============================================================================
// 09vNear: Independent-validation INSPECTOR (template) — NEAR-ONLY variant
// =============================================================================
// =============================================================================
//
// Identical workflow to 09v_valInspector_TEMPLATE.js — only difference is
// FC_ASSET, which points to the 2-strata (S1+S2, no S3) point set produced by
// 09vNear_valGenerator.js + 09vNear_valRandomize.R.
//
// TEMPLATE — each validator copies this to `09vNear_valInspector_<name>.js` so
// the three of them work in parallel and independently.
//
// Workflow per row of your Google-Sheet tab:
//   1. Set PXID to the value in the next row, set YEAR to the year you label.
//   2. Run. Inspect the 1 km neighbourhood: the 4 mosaic layers, the NDYI
//      time-series (Sentinel + HLS), and (in Google Earth Pro) the matching
//      <pxid>.kml for high-resolution historical imagery.
//   3. Fill the row: es_borde, class_2017, class_2017_confidence,
//      class_2025, class_2025_confidence. Switch YEAR to do the other year.
//
// REFERENCE PRIORITY: high-resolution imagery (Google Earth Pro, visible yellow
// bloom) is the PRIMARY evidence. The NDYI/HLS series is phenological SUPPORT
// only — it reuses the model's own features, so don't let it be the sole
// criterion (avoid circularity).
// =============================================================================


// =============================================================================
// 0. VALIDATOR INPUTS  (edit these)
// =============================================================================
var PXID      = 'px001';
var YEAR      = 2017;
var VALIDATOR = 'lican';

var ASSET_PREFIX = 'projects/ee-licanemartinez/assets/Retamap/';
var FC_ASSET     = ASSET_PREFIX + '09vNear_valPoints_fc';

// =============================================================================
// 1. RETRIEVE THE POINT + 1 km NEIGHBOURHOOD
// =============================================================================
var fc      = ee.FeatureCollection(FC_ASSET);
var feat    = ee.Feature(fc.filter(ee.Filter.eq('pxid', PXID)).first());
var contour = feat.geometry();
var centro  = contour.centroid(1);
var aoi     = centro.buffer(500).bounds();   // ~1 km square (clip / region)

Map.centerObject(contour, 18);
Map.setOptions('SATELLITE');
Map.addLayer(aoi,     {color: 'ffffff'}, '1km AOI', false);

// =============================================================================
// 2. FOUR MOSAIC LAYERS from 01_MergedBands_<YEAR>  (NovDec + Feb)
// =============================================================================
var merged = ee.Image(ASSET_PREFIX + '01_MergedBands_' + YEAR).clip(aoi);

var TC_VIS   = {min: 0, max: 3000};
var NDYI_VIS = {min: 0, max: 0.4, palette: ['ffffff', 'ffff66', 'ff9900', 'cc3300']};


var visTC = {
  bands: ['B4', 'B3', 'B2'],
  min: [300, 300, 500],      // Aumenta el 3er valor para reducir el azulado (corta sombras/bruma)
  max: 1400,     // Reduce este valor para aumentar la saturación general
  gamma: 0.6     // Ajusta la gamma si la imagen se ve muy oscura tras el cambio
};

var visTC_feb = Object.assign({}, visTC, {
  bands: ['B4_feb', 'B3_feb', 'B2_feb']
});

Map.addLayer(merged, visTC, YEAR + ' TrueColor NovDec', true);
Map.addLayer(merged.select('NDYI'), NDYI_VIS, YEAR + ' NDYI NovDec', false);
Map.addLayer(merged, visTC_feb, YEAR + ' TrueColor Feb', false);
Map.addLayer(merged.select('NDYI_feb'), NDYI_VIS, YEAR + ' NDYI Feb', false);
Map.addLayer(contour, {color: 'FF0000', fillColor: '00000000'}, 'pixel contour ' + PXID, true);
Map.addLayer(centro,  {color: '00FFFF'}, 'centroid', true);

// =============================================================================
// 3. NDYI TIME SERIES (Aug–Dec of YEAR), clipped to the 1 km AOI
// =============================================================================
var START = ee.Date.fromYMD(YEAR, 8, 1);
var END   = ee.Date.fromYMD(YEAR+1, 3, 31);

// ── 3a. Sentinel-2 (S2_HARMONIZED → no 2022 offset jump) + Cloud Score+ ──────
var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
var s2 = ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
  .filterBounds(aoi)
  .filterDate(START, END)
  .map(function(img) { return img.linkCollection(csPlus, ['cs']); })
  .map(function(img) { return img.updateMask(img.select('cs').gte(0.6)); })
  .map(function(img) {
    return img.normalizedDifference(['B3', 'B2']).rename('NDYI')
      .copyProperties(img, ['system:time_start']);
  });

var chartS2 = ui.Chart.image.series({
  imageCollection: s2.select('NDYI'),
  region: aoi, reducer: ee.Reducer.mean(), scale: 10,
  xProperty: 'system:time_start'
}).setOptions({
  title: 'Sentinel-2 NDYI  Aug–Dec ' + YEAR + '  (' + PXID + ')',
  hAxis: {title: 'date'}, vAxis: {title: 'NDYI (AOI mean)'},
  pointSize: 4, lineWidth: 1
});
// print(chartS2);

// ── 3b. HLS (Landsat + Sentinel harmonised); NDYI is a ratio → no scaling ─────
// Simple Fmask cloud/shadow mask (bits 1=cloud, 3=cloud shadow).
var hlsMask = function(img) {
  var f = img.select('Fmask');
  var clear = f.bitwiseAnd(1 << 1).eq(0)
    .and(f.bitwiseAnd(1 << 3).eq(0));
  return img.updateMask(clear);
};
var hlsNDYI = function(img) {
  return hlsMask(img).normalizedDifference(['B3', 'B2']).rename('NDYI')
    .copyProperties(img, ['system:time_start']);
};
var hlsL = ee.ImageCollection('NASA/HLS/HLSL30/v002')
  .filterBounds(aoi).filterDate(START, END).select(['B2', 'B3', 'Fmask']).map(hlsNDYI);
var hlsS = ee.ImageCollection('NASA/HLS/HLSS30/v002')
  .filterBounds(aoi).filterDate(START, END).select(['B2', 'B3', 'Fmask']).map(hlsNDYI);
var hls = hlsL.merge(hlsS).sort('system:time_start');

var chartHLS = ui.Chart.image.series({
  imageCollection: hls.select('NDYI'),
  region: aoi, reducer: ee.Reducer.mean(), scale: 30,
  xProperty: 'system:time_start'
}).setOptions({
  title: 'HLS (L30+S30) NDYI  Aug–Dec ' + YEAR + '  (' + PXID + ')',
  hAxis: {title: 'date'}, vAxis: {title: 'NDYI (AOI mean)'},
  pointSize: 4, lineWidth: 1, colors: ['1b7837']
});
// print(chartHLS);

// =============================================================================
// 4. ON-SCREEN INSTRUCTIONS
// =============================================================================
print('── VALIDATION ' + PXID + ' · year ' + YEAR + ' · ' + VALIDATOR + ' ──');
print('Primary evidence: high-res imagery in Google Earth Pro (open ' + PXID +
      '.kml). Yellow bloom (Nov–Dec) = retama. NDYI/HLS series = support only.');
print('Fill in your sheet row: es_borde, class_2017, class_2017_confidence, ' +
      'class_2025, class_2025_confidence.');
