// =============================================================================
// 02: RF1 training, prediction, and export to asset
// Runs after 01_generate_sentinel_mosaic.js assets are exported
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER = 'Retamap/Retamap_GEE_Exports';

var exportYears = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];

var roi          = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');
var gt_polygons_new = ee.FeatureCollection(ASSET_PREFIX + 'gt_polys_6_ctrlReduce_moreCtrls_refineRetamas');  

// Band schema from phase 1:
// B2, B3, B4, B8, NDYI, B2_feb, B3_feb, B4_feb, B8_feb, NDYI_feb
var BANDS = ['B2', 'B3', 'B4', 'B8', 'NDYI', 'B2_feb', 'B3_feb', 'B4_feb', 'B8_feb', 'NDYI_feb'];

// =============================================================================
// 1. LOAD mergedCollection FROM ASSETS
// =============================================================================
var mergedCollection = ee.ImageCollection.fromImages(
  exportYears.map(function(year) {
    return ee.Image(ASSET_PREFIX + 'MergedBands_' + year).set('year', year);
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
// 3. STRATIFIED SAMPLING (per year, matched to polygon year)
// =============================================================================
var samples_rf1 = mergedCollection.map(function(image) {
  var year = image.get('year');

  var labeledPolygons = gt_polygons_new
    .filter(ee.Filter.eq('ANIO_FL', year))
    .map(addType01);

  var classImage = labeledPolygons
    .reduceToImage({
      properties: ['type_01'],
      reducer: ee.Reducer.first()
    })
    .rename('type_01');

  return image.addBands(classImage)
    .stratifiedSample({
      numPoints: 5000,
      classBand: 'type_01',
      region: labeledPolygons,
      scale: 10,
      classValues: [0, 1],
      classPoints: [15000, 10000],
      geometries: false,
      tileScale: 16
    })
    .map(function(f) {
      return f.set('year', year);
    });

}).flatten();

// =============================================================================
// 4. TRAIN RF1
// =============================================================================
var rf1 = ee.Classifier.smileRandomForest({numberOfTrees: 500}).train({
  features: samples_rf1,
  classProperty: 'type_01',
  inputProperties: BANDS
});

// =============================================================================
// 5. PREDICT + CONNECTED PIXEL FILTER
// =============================================================================
var predictRF1 = function(image) {
  var pred = image.select(BANDS).classify(rf1).rename('pred');

  // Connected pixel filter applied only to class-1 patches.
  // selfMask() makes class-0 pixels nodata temporarily so connectedPixelCount
  // counts only retama (1) pixels. connected is nodata where pred=0.
  var connected = pred.selfMask().connectedPixelCount(50, true);  // 8-connected

  // Small class-1 patches (< 5 pixels) are reclassified to 0 instead of
  // becoming nodata. Class-0 pixels are preserved as-is.
  // connected.unmask(0): for class-0 pixels (where connected=nodata), unmask
  // gives 0, but pred.eq(1) is false there, so the where() condition is false
  // and class-0 pixels are untouched.
  var predClean = pred
    .where(
      pred.eq(1).and(connected.unmask(0).lt(5)),
      0
    )
    .rename('pred');

  return predClean.set('year', image.get('year'));
};

var preds_collection = mergedCollection.map(predictRF1);

// =============================================================================
// 6. EXPORT RF1 PREDICTIONS TO ASSET
// =============================================================================
exportYears.forEach(function(year) {
  var img = ee.Image(
    preds_collection.filter(ee.Filter.eq('year', year)).first()
  ).select('pred').clip(roi).uint8();

  Export.image.toAsset({
    image: img,
    description: 'Export_RF1_' + year,
    assetId: ASSET_PREFIX + 'RF1_prediction_' + year,
    region: roi,
    scale: 10,
    maxPixels: 1e13
  });

  Export.image.toDrive({
    image         : img,
    description   : 'Drive_RF1_' + year,
    folder        : DRIVE_FOLDER,
    fileNamePrefix: 'RF1_prediction_' + year,
    region        : roi,
    scale         : 10,
    maxPixels     : 1e13,
    fileFormat    : 'GeoTIFF'
  });
});