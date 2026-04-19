// =============================================================================
// 02: RF1 training, prediction, and export to asset
// Runs after 01_sentinelMosaic.js assets are exported
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER    = 'Retamap/Retamap_GEE_Exports';
var EXPORT_TO_DRIVE = false; 
var APPLY_PIXEL_FILTER = true; 
var exportYears = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];

var roi             = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');
var gt_polygons_new = ee.FeatureCollection(ASSET_PREFIX + 'gt_polys_6_ctrlReduce_moreCtrls_refineRetamas');
var compBackgroundPolys = ee.FeatureCollection(ASSET_PREFIX + '04_compBackgroundPolys'); // polys complementarios ctrl
var INCLUDE_COMP_POLYS = true; // Toggle para incluir compBackgroundPolys

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
// 2. LABEL FUNCTION AND MERGE
// =============================================================================
var addType01 = function(feature) {
  var label = ee.Algorithms.If(
    ee.String(feature.get('TIPO')).equals('retama'), 1, 0
  );
  return feature.set('type_01', label);
};

// Formatear polígonos complementarios para asegurar compatibilidad de propiedades
var formatCompPolys = function(feature) {
  // Asignamos type_01 = 0 y un ANIO_FL = 0 (valor dummy, no se filtrará luego)
  return feature.set('type_01', 0).set('ANIO_FL', 0);
};

var labeledGtPolygons = gt_polygons_new.map(addType01);
var labeledPolygons; // Declaración de la variable para el scope global

// Condición controlada por el toggle para incluir o no los polígonos complementarios
if (INCLUDE_COMP_POLYS) {
  var labeledCompPolys = compBackgroundPolys.map(formatCompPolys);
  labeledPolygons = labeledGtPolygons.merge(labeledCompPolys);
} else {
  labeledPolygons = labeledGtPolygons;
}

// =============================================================================
// 3. STRATIFIED SAMPLING
// =============================================================================
var N_CTRL   = 15000;
var N_RETAMA = 10000;

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
  
  // Lógica de filtrado:
  // - Clase 1: extrae solo si ANIO_FL coincide con el año de la imagen.
  // - Clase 0: extrae en todos los años.
  var filterLogic = ee.Filter.or(
    ee.Filter.and(ee.Filter.eq('type_01', 1), ee.Filter.eq('ANIO_FL', year)),
    ee.Filter.eq('type_01', 0)
  );
  
  var yearPoints = samplePoints.filter(filterLogic);

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
    description: 'Export_02_RF1.2_raw_' + year,
    assetId    : ASSET_PREFIX + '02_RF1_2compCTRL_prediction_connectedFilter_' + year,
    region     : roi,
    scale      : 10,
    maxPixels  : 1e13
  });

  if (EXPORT_TO_DRIVE) {
    Export.image.toDrive({
      image         : img,
      description   : 'Drive_02_RF1_raw_' + year,
      folder        : DRIVE_FOLDER,
      fileNamePrefix: '02_RF1_2compCTRL_prediction_connectedFilter_' + year,
      region        : roi,
      scale         : 10,
      maxPixels     : 1e13,
      fileFormat    : 'GeoTIFF'
    });
  }
});