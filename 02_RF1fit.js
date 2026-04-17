// =============================================================================
// 02: RF1 training, prediction, and export to asset
// Runs after 01_sentinelMosaic.js assets are exported
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER    = 'Retamap/Retamap_GEE_Exports';
var EXPORT_TO_DRIVE = false;  // toggle: true → also export to Google Drive

// Nuevo toggle para el filtro de píxeles conectados
var APPLY_PIXEL_FILTER = true; // toggle: true → aplica el filtro para remover parches < 5 px

var exportYears = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];

var roi             = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');
var gt_polygons_new = ee.FeatureCollection(ASSET_PREFIX + 'gt_polys_6_ctrlReduce_moreCtrls_refineRetamas');

// Band schema from 01_sentinelMosaic.js:
var BANDS = ['B2', 'B3', 'B4', 'B8', 'NDYI', 'B2_feb', 'B3_feb', 'B4_feb', 'B8_feb', 'NDYI_feb'];

// =============================================================================
// 1. LOAD mergedCollection FROM ASSETS (01_MergedBands_YYYY)
// =============================================================================
var mergedCollection = ee.ImageCollection.fromImages(
  exportYears.map(function(year) {
    return ee.Image(ASSET_PREFIX + '01_MergedBands_' + year).set('year', year);
  })
);

// =============================================================================
// 2. LABEL FUNCTION
// =============================================================================
var addType01 = function(feature) {
  var label = ee.Algorithms.If(
    ee.String(feature.get('TIPO')).equals('retama'), 1, 0
  );
  return feature.set('type_01', label);
};

// =============================================================================
// 3. STRATIFIED SAMPLING (two-step: fixed N points → per-year pixel lookup)
// =============================================================================
var N_CTRL   = 15000;
var N_RETAMA = 10000;

var labeledPolygons = gt_polygons_new.map(addType01);

var imgType = labeledPolygons.reduceToImage({
  properties: ['type_01'],
  reducer: ee.Reducer.first()
}).rename('type_01');

var imgAnioFL = labeledPolygons.reduceToImage({
  properties: ['ANIO_FL'],
  reducer: ee.Reducer.first()
}).rename('ANIO_FL');

var labelImage = imgType.addBands(imgAnioFL);

var samplePoints = labelImage
  .stratifiedSample({
    numPoints  : 0,
    classBand  : 'type_01',
    region     : roi,
    scale      : 10,
    classValues: [0, 1],
    classPoints: [N_CTRL, N_RETAMA],
    geometries : true,
    tileScale  : 16
  });

var samples_rf1 = mergedCollection.map(function(image) {
  var year = image.get('year');
  var yearPoints = samplePoints.filter(ee.Filter.eq('ANIO_FL', year));

  return image.select(BANDS)
    .sampleRegions({
      collection : yearPoints,
      properties : ['type_01', 'ANIO_FL'],
      scale      : 10,
      geometries : false
    })
    .map(function(f) { return f.set('year', year); });
}).flatten();

// =============================================================================
// 4. TRAIN RF1
// =============================================================================
var rf1 = ee.Classifier.smileRandomForest({numberOfTrees: 500}).train({
  features       : samples_rf1,
  classProperty  : 'type_01',
  inputProperties: BANDS
});

// =============================================================================
// 5. PREDICT (+ OPTIONAL CONNECTED PIXEL FILTER)
// =============================================================================
var predictRF1 = function(image) {
  var pred = image.select(BANDS).classify(rf1).rename('pred');

  // Evalúa el toggle cliente para decidir si aplicar el filtrado o no
  if (APPLY_PIXEL_FILTER) {
    var connected = pred.selfMask().connectedPixelCount(50, true);
    pred = pred.where(
      pred.eq(1).and(connected.unmask(0).lt(5)),
      0
    ).rename('pred');
  }

  return pred.set('year', image.get('year'));
};

var preds_collection = mergedCollection.map(predictRF1);

// =============================================================================
// 6. EXPORT RF1 PREDICTIONS TO ASSET (02_RF1_raw_prediction_YYYY)
// =============================================================================
exportYears.forEach(function(year) {
  var img = ee.Image(
    preds_collection.filter(ee.Filter.eq('year', year)).first()
  ).select('pred').clip(roi).uint8();

  Export.image.toAsset({
    image      : img,
    description: 'Export_02_RF1_raw_' + year,
    assetId    : ASSET_PREFIX + '02_RF1_prediction_connectedFilter_' + year,
    region     : roi,
    scale      : 10,
    maxPixels  : 1e13
  });

  if (EXPORT_TO_DRIVE) {
    Export.image.toDrive({
      image         : img,
      description   : 'Drive_02_RF1_raw_' + year,
      folder        : DRIVE_FOLDER,
      fileNamePrefix: '02_RF1_prediction_connectedFilter_' + year,
      region        : roi,
      scale         : 10,
      maxPixels     : 1e13,
      fileFormat    : 'GeoTIFF'
    });
  }
});