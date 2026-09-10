

// 9-Invasion_area_calculation.js
// Calculo de la diferencia en cobertura de invasion de Retama
// -------------------------------

// -------------------------------
// Geometrías y datos base

// -------------------------------
var studyArea = table("projects/ee-sofiapptemp1/assets/3_study_area_retama")
var roi = studyArea

Map.addLayer(roi, {palette:['black']}, 'ROI', false);

// Infraestructura y capas auxiliares
var roads = ee.FeatureCollection("users/IvanBarbera/patagonian_fires/roads");
var humans = ee.FeatureCollection("projects/cosmic-anthem-398720/assets/3_poblados_study_area");
var lines_water = ee.FeatureCollection("projects/cosmic-anthem-398720/assets/lineas_de_aguas_continentales_perenne");
var areas_water = ee.FeatureCollection("projects/cosmic-anthem-398720/assets/areas_de_aguas_continentales_perenne");
var apn = ee.FeatureCollection("projects/cosmic-anthem-398720/assets/apn_selected_PatagoniaNorte");

// Vegetación y elevación
var vegetation = ee.Image("users/IvanBarbera/Fire_spread/vegetation_ciefap_wwf");
var nasa_elevation = ee.Image("NASA/NASADEM_HGT/001");


//RF results
var years = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];

// Creamos una lista de imágenes llamando a cada asset por su ruta individual
var imageList = years.map(function(year) {
  var assetPath = 'projects/cosmic-anthem-398720/assets/RF_retama_raster_years/08_RF2_prediction_' + year; // Ajusta el nombre exacto del asset
  return ee.Image(assetPath).set('year', year);
});

// Convertimos esa lista en una ImageCollection real
var RFresults = ee.ImageCollection.fromImages(imageList);
var YearForVis = 2023


// Rango temporal
var startYear = 2017;
var endYear = 2025;


// Creacion de capas para hacer mediciones de area

// --- Buffer de 20m para rutas (Eje central + 20m a cada lado) ---
var roadsBuffer = roads.distance(20).lte(20);

// --- Buffer de 20m para rios (Eje central + 20m a cada lado) ---

var waterLinesBuffer = lines_water.distance(20).lte(20);

// --- Buffer de "Costa" para Cuerpos de Agua (Solo los 20m exteriores) ---
// 1. Convertimos los polígonos de agua en una imagen binaria
var waterAreasImg = ee.Image().paint(areas_water, 1);

// 2. Creamos un buffer de 20m alrededor de los polígonos
var waterBufferAll = areas_water.distance(20).lte(20);

// 3. RESTAMOS el agua original para quedarnos solo con la "tierra" de la orilla
// .and(waterAreasImg.not()) elimina el interior del lago
var coastalBuffer = waterBufferAll.and(waterAreasImg.not());
var coastalBufferMask = coastalBuffer.selfMask();

// --- Urbanizaciones y No Urbanizaciones (limitado al ROI) ---
var urbanMask = ee.Image(0).paint(humans, 1).selfMask().clip(roi);
var nonUrbanMask = urbanMask.unmask(0).not().selfMask().clip(roi);

// --- APN y Fuera de APN (dentro del ROI) ---
var apnMask = ee.Image().paint(apn, 1); 
var nonApnMask = apnMask.unmask(0).not().selfMask().clip(roi);

//funcion para calcular Hectarias

var calculateArea = function(classification, mask, label) {
  // 1. Aseguramos que trabajamos con una imagen de 1 banda y aplicamos la máscara
  var areaImage = classification
    .select(0) // Selecciona la primera banda sea cual sea el nombre
    .selfMask() 
    .updateMask(mask)
    .multiply(ee.Image.pixelArea()); // Área en m²

  var stats = areaImage.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: roi,
    scale: 10, 
    maxPixels: 1e13
  });

// 2. Extraemos el valor usando la primera clave disponible en el diccionario
  var areaM2 = ee.Number(stats.get(stats.keys().get(0)));
  
  // 3. Manejo de casos donde el área sea 0 (evita que devuelva null)
  var areaHa = ee.Algorithms.If(areaM2, areaM2.divide(10000), 0);

  return ee.Feature(null, {
    'zona': label,
    'hectareas': areaHa,
    'anio': 2023
  });
};

// Seleccionamos la imagen de 2023 de tu colección
var classification2023 = RFresults.filter(ee.Filter.eq('year', 2023)).first();


// Lista de cálculos
var areaResults = ee.FeatureCollection([
  calculateArea(classification2023, urbanMask, "Urbanizaciones"),
  calculateArea(classification2023, nonUrbanMask, "Fuera de Urbanizaciones"),
  calculateArea(classification2023, roadsBuffer, "Buffer Rutas 20m"),
  calculateArea(classification2023, waterLinesBuffer, "Buffer rios 20m"),
  calculateArea(classification2023, coastalBuffer, "Buffer costas de lagos"),
  calculateArea(classification2023, apnMask, "Dentro APN"),
  calculateArea(classification2023, nonApnMask, "Fuera APN")
]);

print("Resultados de Área 2023 (Ha):", areaResults);

// Exportar a CSV si lo necesitas
Export.table.toDrive({
  collection: areaResults,
  description: 'Areas_Retama_2023_Subzonas',
  fileFormat: 'CSV'
});

// Visualizar para chequear solapamiento
Map.addLayer(urbanMask, {palette:['red']}, 'Máscara Urbana', false);

Map.addLayer(roadsBuffer, {palette:['green']}, 'Máscara rutas', false);
Map.addLayer(classification2023.selfMask(), {palette:['yellow']}, 'Retama 2023', false);

// ================================================================
// CÁLCULO DE % DE OCUPACIÓN para 1 año
// ================================================================

// 1. Definimos la función de cálculo una sola vez
var calculateOccupancy = function(classification, mask, label) {
  // Área total de la zona (ej: toda la superficie de rutas)
  var totalArea = ee.Image.pixelArea().updateMask(mask).reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: roi,
    scale: 10,
    maxPixels: 1e13
  });
  var totalHa = ee.Number(totalArea.get('area', 0)).divide(10000);

  // Área de Retama detectada dentro de esa zona
  var retamaArea = ee.Image.pixelArea().updateMask(mask).updateMask(classification.gt(0)).reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: roi,
    scale: 10,
    maxPixels: 1e13
  });
  var retamaHa = ee.Number(retamaArea.get('area', 0)).divide(10000);

  // Porcentaje: (Retama / Total Zona) * 100
  var porcentaje = ee.Number(ee.Algorithms.If(totalHa.gt(0), 
    retamaHa.divide(totalHa).multiply(100), 0));

  return ee.Feature(null, {
    'zona': label,
    'ha_total_zona': totalHa,
    'ha_retama': retamaHa,
    'porcentaje_ocupacion': porcentaje,
    'anio': 2023
  });
};

// 2. Ejecutamos todos los cálculos en una sola lista (más ordenado)
var occupancyResults2023 = ee.FeatureCollection([
  calculateOccupancy(classification2023, urbanMask, "Urbanizaciones"),
  calculateOccupancy(classification2023, nonUrbanMask, "Fuera de Urbanizaciones"),
  calculateOccupancy(classification2023, roadsBuffer, "Buffer Rutas 20m"),
  calculateOccupancy(classification2023, waterLinesBuffer, "Buffer rios 20m"),
  calculateOccupancy(classification2023, coastalBufferMask, "Buffer costas de lagos"),
  calculateOccupancy(classification2023, apnMask, "Dentro APN"),
  calculateOccupancy(classification2023, nonApnMask, "Fuera APN")
]);

print("Reporte Final Ocupación 2023:", occupancyResults2023);

// Exportar este reporte específico si quieres
Export.table.toDrive({
  collection: occupancyResults2023,
  description: 'Porcentaje_Ocupacion_Retama_2023',
  fileFormat: 'CSV'
});

// ================================================================
// CÁLCULO DE % DE OCUPACIÓN para todos los años
// ================================================================
var zones = ee.Image(0)
  .where(urbanMask, 1)
  .where(nonUrbanMask, 2)
  .where(roadsBuffer, 3)
  .where(waterLinesBuffer, 4)
  .where(coastalBufferMask, 5)
  .where(apnMask, 6)
  .where(nonApnMask, 7)
  .rename('zona');

Map.addLayer(zones.randomVisualizer(), {}, "Zonas", false);
var calculateOccupancyFast = function(image){

  var year = image.get('year');

  var retama = image.select(0).gt(0);

  var areaImage = ee.Image.pixelArea().divide(10000);

  var totalArea = areaImage.addBands(zones)
    .reduceRegion({
      reducer: ee.Reducer.sum().group({
        groupField: 1,
        groupName: 'zona'
      }),
      geometry: roi,
      scale: 10,
      maxPixels: 1e13
    });

  var retamaArea = areaImage.updateMask(retama)
    .addBands(zones)
    .reduceRegion({
      reducer: ee.Reducer.sum().group({
        groupField: 1,
        groupName: 'zona'
      }),
      geometry: roi,
      scale: 10,
      maxPixels: 1e13
    });

  var totalDict = ee.List(totalArea.get('groups'));
  var retamaDict = ee.List(retamaArea.get('groups'));

  var results = totalDict.map(function(item){

    item = ee.Dictionary(item);

    var zona = item.get('zona');
    var totalHa = ee.Number(item.get('sum'));

    var retItem = ee.Dictionary(
      retamaDict.filter(ee.Filter.eq('zona', zona)).get(0)
    );

    var retHa = ee.Number(retItem.get('sum', 0));

    var perc = ee.Algorithms.If(
      totalHa.gt(0),
      retHa.divide(totalHa).multiply(100),
      0
    );

    return ee.Feature(null,{
      zona: zona,
      ha_total_zona: totalHa,
      ha_retama: retHa,
      porcentaje_ocupacion: perc,
      anio: year
    });

  });

  return ee.FeatureCollection(results);
};

var occupancyAllYears = RFresults
  .map(calculateOccupancyFast)
  .flatten();

// print("Ocupación todos los años:", occupancyAllYears); 
// El export anda, esto se sobrepasa de tiempo de calculo

Export.table.toDrive({
  collection: occupancyAllYears,
  description: 'Porcentaje_Ocupacion_Retama_2017_2023',
  fileFormat: 'CSV'
});

// ================================================================
// CÁLCULO DE ÁREA TOTAL DE RETAMA (Todo el ROI) 2017-2023
// ================================================================
var calculateTotalArea = function(image) {
  var year = image.get('year');
  
  // 1. Forzamos a que CUALQUIER píxel con retama valga 1
  // gt(0) convierte todo lo que no es cero en 1.
  var retamaBinaria = image.select(0).gt(0).selfMask();
  
  // 2. Usamos el reductor COUNT para saber cuántos píxeles hay exactamente
  var stats = retamaBinaria.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: roi,
    scale: 10, 
    maxPixels: 1e13
  });

  var count = ee.Number(stats.get(stats.keys().get(0)));
  
  // 3. Calculamos el área: (Cantidad de píxeles * área de 1 píxel) / 10000
  // El área de un píxel a escala de sentinel es 75.14m²
  var areaHa = count.multiply(75.14).divide(10000); 

  return ee.Feature(null, {
    'anio': year,
    'area_total_ha': areaHa,
    'tipo': 'Total ROI'
  });
};

// Aplicar a todos los años
var totalAreaHistory = RFresults.map(calculateTotalArea);

// 1. Mostrar resultados en la Consola
print("Evolución Área Total Retama (Ha):", totalAreaHistory);

// 2. Crear Gráfico con .groups (La forma más estable)
var totalChart = ui.Chart.feature.groups({
  features: totalAreaHistory,
  xProperty: 'anio',
  yProperty: 'area_total_ha',
  seriesProperty: 'tipo'
}).setChartType('ColumnChart') // O 'LineChart' si prefieres líneas
  .setOptions({
    title: 'Superficie Total de Retama por Año (Hectáreas)',
    hAxis: {title: 'Año', format: '####'},
    vAxis: {title: 'Hectáreas'},
    colors: ['#e34a33']
  });

print(totalChart);

// 3. Exportar a CSV
Export.table.toDrive({
  collection: totalAreaHistory,
  description: 'Retama_Area_Total_2017_2023',
  fileFormat: 'CSV'
});

// ================================================================
// REFERENCIA VISUAL: DIBUJAR 1 HECTÁREA EN EL MAPA
// ================================================================

// 1. Definimos un punto central para la referencia 
// (Puedes cambiar estas coordenadas por unas de tu zona de interés)
var puntoReferencia = ee.Geometry.Point([-71.41232742499766,  -41.12133799985396]); // Cerca de Bariloche

// 2. Creamos un buffer cuadrado de 50m de radio (lo que genera un cuadrado de 100x100m)
var unaHectarea = puntoReferencia.buffer(50).bounds();

// 3. Agregamos la capa al mapa con un color llamativo (azul)
Map.addLayer(unaHectarea, {color: '00FFFF'}, 'REFERENCIA: 1 Hectárea (100x100m)');

// 4. Centramos el mapa en esa hectárea para verla de cerca
Map.centerObject(puntoReferencia, 16);

print("Referencia: El cuadrado azul en el mapa representa exactamente 1 Hectárea.");

// ================================================================
// TRANSICIONES de Area ANUALES
// ================================================================

// --- Función para calcular transiciones entre dos años ---
var getTransition = function(year1, year2) {
  var img1 = RFresults.filter(ee.Filter.eq('year', year1)).first().select(0).gt(0);
  var img2 = RFresults.filter(ee.Filter.eq('year', year2)).first().select(0).gt(0);
  
  // Lógica de transición:
  // Si habia en ambos = 1 (Permanencia - Gris)
  // Si no habia en 1 y hay en 2 = 2 (Ganancia - Rojo)
  // Si habia en 1 y no hay en 2 = 3 (Pérdida - Azul)
  
  var transition = ee.Image(0)
    .where(img1.and(img2), 1)
    .where(img1.not().and(img2), 2)
    .where(img1.and(img2.not()), 3)
    .selfMask() // Quitamos los ceros (donde nunca hubo nada)
    .clip(roi);
    
  return transition.set('transition_period', year1 + '-' + year2);
};

// --- Generamos las 6 transiciones consecutivas ---
var trans17_18 = getTransition(2017, 2018);
var trans18_19 = getTransition(2018, 2019);
var trans19_20 = getTransition(2019, 2020);
var trans20_21 = getTransition(2020, 2021);
var trans21_22 = getTransition(2021, 2022);
var trans22_23 = getTransition(2022, 2023);
var trans22_23 = getTransition(2023, 2024);
var trans22_23 = getTransition(2024, 2025);
// -- genero una transicion general
var trans17_25 = getTransition(2017, 2025);


// vizualizar
var transVis = {
  min: 1,
  max: 3,
  palette: [
    '#F2F1E4', // 1: Gris (Estable)
    '#E64E4B', // 2: Rojo (Invasión nueva)
    '#023FA5'  // 3: Azul (Desapareció/Pérdida)
  ]
};

// Agregamos, por ejemplo, los últimos dos periodos
Map.addLayer(trans17_18, transVis, 'Transición 2017-2018');
Map.addLayer(trans19_20, transVis, 'Transición 2019-2020');
Map.addLayer(trans22_23, transVis, 'Transición 2022-2023');

Map.addLayer(trans17_25, transVis, 'Transición 2017-2025');


// Función para calcular áreas de transición corregida
var calculateTransitionArea = function(year1, year2) {
  var period = ee.String(ee.Number(year1).format()).cat('-').cat(ee.Number(year2).format());
  var transImage = getTransition(year1, year2);
  
  var stats = transImage.reduceRegion({
    reducer: ee.Reducer.frequencyHistogram(),
    geometry: roi,
    scale: 10,
    maxPixels: 1e13
  });

  // Obtenemos el diccionario (ej: {"1": 450, "2": 120...})
  var counts = ee.Dictionary(stats.get(stats.keys().get(0)));
  
  // Función interna para generar los Features de cada clase
  var labels = [['1', 'Estable'], ['2', 'Ganancia (Invasión)'], ['3', 'Pérdida (Control)']];
  
  var features = labels.map(function(item) {
    var val = ee.String(item[0]);
    var label = ee.String(item[1]);
    
    // Obtenemos el conteo del diccionario. Si no existe la clave, devolvemos 0.
    var count = ee.Number(counts.get(val, 0)); 
    
    var areaHa = count.multiply(100).divide(10000);
    
    return ee.Feature(null, {
      'periodo': period,
      'tipo': label,
      'hectareas': areaHa
    });
  });

  return ee.FeatureCollection(features);
};

// --- Generamos todas las transiciones automáticamente ---


// --- Gráfico de Transiciones ---

// --- 1. Ajuste de la función para asegurar que 'periodo' no se pierda ---
var calculateTransitionStats = function(year1, year2) {
  var y1 = ee.Number(year1);
  var y2 = ee.Number(year2);
  var period = ee.String(y1.format()).cat('-').cat(y2.format());
  
  var transImg = getTransition(y1, y2);
  
  var stats = transImg.reduceRegion({
    reducer: ee.Reducer.frequencyHistogram(),
    geometry: roi,
    scale: 10,
    maxPixels: 1e13
  });

  var counts = ee.Dictionary(stats.get(stats.keys().get(0)));
  
  // Extraemos conteos asegurando que si no existen sean 0
  var c1 = ee.Number(counts.get('1', 0)); // Estable
  var c2 = ee.Number(counts.get('2', 0)); // Ganancia
  var c3 = ee.Number(counts.get('3', 0)); // Pérdida

  // Creamos los 3 Features manualmente para asegurar que tengan la propiedad 'periodo'
  var f1 = ee.Feature(null, {'periodo': period, 'tipo': 'Estable', 'hectareas': c1.multiply(100).divide(10000)});
  var f2 = ee.Feature(null, {'periodo': period, 'tipo': 'Ganancia (Invasión)', 'hectareas': c2.multiply(100).divide(10000)});
  var f3 = ee.Feature(null, {'periodo': period, 'tipo': 'Pérdida (Control)', 'hectareas': c3.multiply(100).divide(10000)});

  return ee.FeatureCollection([f1, f2, f3]);
};

// --- 2. Generar los datos para todos los años ---
var yearsList = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var transitionStats = ee.FeatureCollection(
  yearsList.slice(0, -1).map(function(y) {
    return calculateTransitionStats(y, y + 1);
  })
).flatten();

// --- 3. Imprimir y Graficar ---
print('Verificación de Tabla (mira si aparece periodo aquí):', transitionStats.limit(5));

var transitionChart = ui.Chart.feature.groups({
  features: transitionStats,
  xProperty: 'periodo',
  yProperty: 'hectareas',
  seriesProperty: 'tipo'
}).setChartType('ColumnChart')
  .setOptions({
    title: 'Dinámica de Cambio Anual: Ganancia vs Pérdida',
    hAxis: {title: 'Periodo de Transición'},
    vAxis: {title: 'Hectáreas'},
    colors: ['#F2F1E4', '#E64E4B', '#023FA5']
  });

print(transitionChart);
