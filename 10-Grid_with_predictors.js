var grid_1000_empty = ee.FeatureCollection("projects/cosmic-anthem-398720/assets/grid_1000m_clipped_vacio"),
    studyArea = ee.FeatureCollection("projects/ee-sofiapptemp1/assets/3_study_area_retama"),
    grid_7000_empty = ee.FeatureCollection("projects/cosmic-anthem-398720/assets/grid_7000m_clipped_vacio"),
    grid_3000_empty = ee.FeatureCollection("projects/cosmic-anthem-398720/assets/grid_3000m_clipped_vacio"),
    years_retama = ee.FeatureCollection("projects/cosmic-anthem-398720/assets/RF_retama_raster_years/retama_all_years"),
    img_years_retama = ee.Image("projects/cosmic-anthem-398720/assets/RF_retama_raster_years/retama_raster_stack");
    
/**
 * Análisis de Predictores en Grilla (Patagonia)
 * --------------------------------------------
 * Este script calcula variables espaciales y climáticas
 * para cada celda de una grilla de 1000 * 1000 m.
 */

// ==============================================================================
// 1. IMPORTACIÓN Y CONFIGURACIÓN DE DATOS
// ==============================================================================

// Máscara de área disponible
var availableArea = ee.Image("projects/cosmic-anthem-398720/assets/available_area_mask").selfMask();

// Capas vectoriales
var roads = ee.FeatureCollection("users/IvanBarbera/patagonian_fires/roads");
var fires = ee.FeatureCollection("projects/cosmic-anthem-398720/assets/patagonian_fires");
var linesWater = ee.FeatureCollection("projects/cosmic-anthem-398720/assets/lineas_de_aguas_continentales_perenne");
var areasWater = ee.FeatureCollection("projects/cosmic-anthem-398720/assets/1_lagos_study_area");
var apn = ee.FeatureCollection("projects/cosmic-anthem-398720/assets/apn_selected_PatagoniaNorte");
var coastLines = ee.FeatureCollection("projects/cosmic-anthem-398720/assets/2_coastLines_study_area");

// Capas raster
var landcovers = ee.Image("projects/mapbiomas-argentina/assets/LAND-COVER/COLLECTION-2/INTEGRATION/mapbiomas_argentina_collection1_integration_v8");
var vegetation = ee.Image("users/IvanBarbera/Fire_spread/vegetation_ciefap_wwf");
var nasaElevation = ee.Image("NASA/NASADEM_HGT/001");

// ==============================================================================
// 2. PRE-PROCESAMIENTO DE CAPAS
// ==============================================================================

// Conversión a imágenes binarias (Masks)
var firesImg = ee.Image().byte().paint(fires, 1).selfMask();
var apnImg = ee.Image().byte().paint(apn, 1).selfMask();
var nonApnMask = apnImg.unmask(0).not().selfMask().clip(studyArea);

var waterAreasImg = ee.Image().byte().paint(areasWater, 1).selfMask();
var waterMask = waterAreasImg.gt(0).selfMask();

// Capas urbano
print('Bandas disponibles en LandCover:', landcovers.bandNames());
// Definir el rango de años
var startYear = 2017;
var endYear = 2024; 
var years = [];
for (var i = startYear; i <= endYear; i++) {
  years.push(i);
}

// urbano 2023
var urbano2023 = landcovers.select('classification_2023').eq(24).rename('urbano').updateMask(availableArea);
var nonUrbanMask = urbano2023.unmask(0).not().selfMask().clip(studyArea);

// Climatología (TerraClimate)
var climate = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filterDate('2017-01-01', '2024-12-31');
var precipMean = climate.select('pr').mean();
var tmean = climate.map(function(img){
  return img.expression('(tmmx + tmmn) / 2', {tmmx: img.select('tmmx'), tmmn: img.select('tmmn')}).rename('tmean');
}).mean().divide(10);
var tmaxMean = climate.select('tmmx').mean().divide(10);
var tminMean = climate.select('tmmn').mean().divide(10);

// retama bandas disponibles
print('Bandas disponibles en retama:', img_years_retama.bandNames());
// Lista de años

// Generar nombres tipo "retama_2017", "retama_2018", ...
var newNames = ee.List.sequence(2017, 2025).map(function(y) {
  return ee.String('retama_').cat(ee.Number(y).format('%d'));
});

var img_renamed = img_years_retama.rename(newNames);

print(img_renamed.bandNames());


var bandNames = img_years_retama.bandNames();

var retamaCollection = ee.ImageCollection(
  bandNames.map(function(bandName) {
    
    var year = ee.Number.parse(
      ee.String(bandName).replace('b', '')
    ).add(2016); // si b1 = 2017
    
    return img_years_retama
      .select([bandName])
      .rename('retama')
      .set({
        'year': year,
        'system:time_start': ee.Date.fromYMD(year, 1, 1).millis()
      });
  })
);
print(retamaCollection);

// ==============================================================================
// 3. FUNCIONES AUXILIARES
// ==============================================================================

/** Calcula porcentaje de cobertura de una máscara en una geometría */
var calcPercentage = function(feature, maskImage) {
  var areaImage = ee.Image.pixelArea();
  var masked = areaImage.updateMask(maskImage).updateMask(availableArea).reduceRegion({
    reducer: ee.Reducer.sum(), 
    geometry: feature.geometry(), 
    scale: 10, 
    maxPixels: 1e13,
    tileScale: 4
  }).get('area');

  var total = areaImage.updateMask(availableArea).reduceRegion({
    reducer: ee.Reducer.sum(), 
    geometry: 
    feature.geometry(), 
    scale: 10, 
    maxPixels: 1e13,
    tileScale: 4
  }).get('area');

  return ee.Algorithms.If(ee.Number(total).gt(0), ee.Number(masked).divide(ee.Number(total)).multiply(100), 0);
};

/** Calcula longitud de líneas dentro de una geometría */
var calcLength = function(feature, linesFC) {
  var intersect = linesFC.filterBounds(feature.geometry()).map(function(f){
    return f.intersection(feature.geometry(), 1);
  });
  return ee.Algorithms.If(intersect.size().gt(0), intersect.geometry().length(), 0);
};


/** Extrae valor medio de una imagen climática */
var extractMean = function(image, bandName, cell) {
  var stats = image.reduceRegion({
    reducer: ee.Reducer.mean(), 
    geometry: cell.geometry(), 
    scale: 4000, 
    maxPixels: 1e13,
    tileScale: 4
  });
  return ee.Number(stats.get(bandName, 0)); // Retorna 0 si es nulo
};

// ==============================================================================
// 4. PROCESAMIENTO DE GRILLA
// ==============================================================================

// Filtrar celdas válidas (> 100 m2)
var gridValid = grid_1000_empty.map(function(cell) {
  var area = ee.Image.pixelArea().updateMask(availableArea).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: cell.geometry(), scale: 10, maxPixels: 1e13
  });
  return cell.set('available_area_m2', ee.Number(area.get('area', 0)));
}).filter(ee.Filter.gt('available_area_m2', 1000));

// Calcular variables para cada celda
var result = gridValid.map(function(cell) {
  var centroid = cell.geometry().centroid(1).coordinates();
  
  // -- Inicio bloque urbano --
  var urbanData = {}; // Objeto vacío para guardar los valores
  
  years.forEach(function(year) {
    var bandName = 'classification_' + year;
    
    // Verificamos si la banda existe para evitar errores
    if (landcovers.bandNames().contains(bandName)) {
      var img = landcovers.select(bandName).eq(24).rename('urbano').updateMask(availableArea);
      urbanData['urbano_perc_' + year] = calcPercentage(cell, img);
    }
  });
  // -- Fin del bloque --

// -- Inicio bloque urbano --
var retamaData = {};

years.forEach(function(year) {
  
  var img = retamaCollection
    .filter(ee.Filter.eq('year', year))
    .first(); // obtenés la imagen de ese año
  
  // chequeo por si falta algún año
  if (img) {
    var retamaMask = img
      .eq(1) // asumimos 1 = presencia de retama
      .rename('retama')
      .updateMask(availableArea);
    
    retamaData['retama_perc_' + year] = calcPercentage(cell, retamaMask);
  }
});
  // -- Fin del bloque --
  
  // Construimos el objeto final
  var properties = {
    lon: centroid.get(0),
    lat: centroid.get(1),
    // ... tus otras variables (incendios, apn, etc.)
    roads_len: calcLength(cell, roads),
    rivers_len: calcLength(cell, linesWater),
    coast_len: calcLength(cell, coastLines),
    apn_perc: calcPercentage(cell, apnImg),
    non_apn_perc: calcPercentage(cell, nonApnMask),
    fires_perc: calcPercentage(cell, firesImg),
    // Clima
    elev_mean: extractMean(nasaElevation, 'elevation', cell),
    precip_mean: extractMean(precipMean, 'pr', cell),
    temp_mean: extractMean(tmean, 'tmean', cell)
    
  };
    
 return cell.set(properties).set(urbanData).set(retamaData);
});
// Vizualizar resultado

// 1. Define el ID que buscas
var idBuscado = '00000000000000003090'; // Cambia esto por el número de grilla que quieres consultar

// 2. Filtra la colección
var celdaEspecifica = result.filter(ee.Filter.eq('system:index', idBuscado));

// 3. Imprime el resultado para ver sus valores
print('Datos de la celda ' + idBuscado + ':', celdaEspecifica);

// ==============================================================================
// 5. EXPORTACIÓN Y VISUALIZACIÓN
// ==============================================================================

/* Export.table.toDrive({
  collection: result,
  description: 'grid_1000_predictors_table',
  fileFormat: 'CSV'
});

Map.centerObject(grid_1000_empty, 12);
Map.addLayer(grid_1000_empty, {color: 'blue'}, 'Mi Grilla de 1000');
Map.addLayer(grid_3000_empty, {color: 'yellow'}, 'Mi Grilla de 3000');
Map.addLayer(grid_7000_empty, {color: 'green'}, 'Mi Grilla de 7000');
*/

// ==============================================================================
//funcion para cualquier grilla
// ==============================================================================

// Esta función procesa CUALQUIER grilla que le pases
var processGrid = function(inputGrid, threshold) {
  
  // 1. Filtrar celdas válidas (usamos el inputGrid recibido)
  var gridValid = inputGrid.map(function(cell) {
    var area = ee.Image.pixelArea().updateMask(availableArea).reduceRegion({
      reducer: ee.Reducer.sum(), geometry: cell.geometry(), scale: 10, maxPixels: 1e13
    });
    return cell.set('available_area_m2', ee.Number(area.get('area', 0)));
  }).filter(ee.Filter.gt('available_area_m2', threshold)); // Ajusta este número si es necesario

  // 2. Calcular variables
  var result = gridValid.map(function(cell) {
    var centroid = cell.geometry().centroid(1).coordinates();
    var urbanData = {};
    var retamaData = {};
    
    years.forEach(function(year) {
      var bandName = 'classification_' + year;
      if (landcovers.bandNames().contains(bandName)) {
        var img = landcovers.select(bandName).eq(24).rename('urbano').updateMask(availableArea);
        var urbPerc = calcPercentage(cell, img);
        urbanData['urbano_perc_' + year] = urbPerc;
        urbanData['no_urbano_perc_' + year] = ee.Number(100).subtract(urbPerc);
      }
    });
    
    
    years.forEach(function(year) {
  
    
  var img = retamaCollection
    .filter(ee.Filter.eq('year', year))
    .first();
  
  // chequeo por si ese año no existe
  if (img) {
    
    var retamaMask = img
      .eq(1) 
      .rename('retama')
      .updateMask(availableArea);
    
    var retPerc = calcPercentage(cell, retamaMask);
    
    retamaData['retama_perc_' + year] = retPerc;
  }
});
    
    
    var properties = {
      lon: centroid.get(0),
      lat: centroid.get(1),
      roads_len: calcLength(cell, roads),
      rivers_len: calcLength(cell, linesWater),
      coast_len: calcLength(cell,coastLines),
      apn_perc: calcPercentage(cell, apnImg),
      non_apn_perc: calcPercentage(cell, nonApnMask),
      fires_perc: calcPercentage(cell, firesImg),
      elev_mean: extractMean(nasaElevation, 'elevation', cell),
      precip_mean: extractMean(precipMean, 'pr', cell),
      temp_mean: extractMean(tmean, 'tmean', cell)
    };
    return cell.set(properties).set(urbanData).set(retamaData);
  });
  
  return result;
};

// Procesa cada una
var result1000 = processGrid(grid_1000_empty, 1000);
var result3000 = processGrid(grid_3000_empty, 1000);
var result7000 = processGrid(grid_7000_empty, 1000);

// Exporta cada una
Export.table.toDrive({ collection: result1000, description: 'grid_1000_predictors', fileFormat: 'CSV' });
Export.table.toDrive({ collection: result3000, description: 'grid_3000_predictors', fileFormat: 'CSV' });
Export.table.toDrive({ collection: result7000, description: 'grid_7000_predictors', fileFormat: 'CSV' });


Map.addLayer(result1000, {color: 'blue'}, 'Mi Grilla de 1000');
Map.addLayer(result3000, {color: 'blue'}, 'Mi Grilla de 3000');
Map.addLayer(result7000, {color: 'blue'}, 'Mi Grilla de 7000');