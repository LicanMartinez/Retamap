// =============================================================================
// 13: Grid invasion proportion — per-cell retama pixel count and proportion
// -----------------------------------------------------------------------------
// Prerequisite : asset 12_grid_{CELL_SIZE}m_clipped  (from 12_gridExport.js)
//              + assets 08_RF2_prediction_YYYY{SUFFIX} (from 08_RF2patchFilter.js)
// Produces     : asset 13_gridInvasion_{CELL_SIZE}m{SUFFIX}
//              + Drive shapefile (if EXPORT_TO_DRIVE = true)
// -----------------------------------------------------------------------------
// MOD 4: For each grid cell and each year, computes:
//   retama_YYYY   — sum of class-1 pixels (= retama pixel count)
//   total_YYYY    — total non-masked pixel count in that cell
//   proportion_YYYY — retama_YYYY / total_YYYY
//
// Note: sum() works because classification is 0/1 → sum = count of class-1.
// count() gives total non-masked pixels, handling border cells and masked pixels
// (water, elevation) correctly — proportion is retama/valid, not retama/area.
// =============================================================================

// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER    = 'Retamap/Retamap_GEE_Exports';
var EXPORT_TO_DRIVE = false;  // toggle: true → also export shapefile to Drive

var exportYears = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025];
var CELL_SIZE   = 3000;   // metres — must match the grid exported by 12_gridExport.js
var RUN_SUFFIX  = '_s1.10k_qs1.0_s0.20k_qs0.0';

// =============================================================================
// 1. LOAD GRID
// =============================================================================
var grid = ee.FeatureCollection(ASSET_PREFIX + '12_grid_' + CELL_SIZE + 'm_clipped');
print('Grid cell count:', grid.size());

// =============================================================================
// 2. ITERATIVELY ADD PER-YEAR STATS VIA reduceRegions
// =============================================================================
var gridResult = grid;

exportYears.forEach(function(year) {
  var pred = ee.Image(ASSET_PREFIX + '08_RF2_prediction_' + year + RUN_SUFFIX)
    .select('classification');

  // sum = count of class-1 pixels (retama); count = total non-masked pixels
  gridResult = pred.reduceRegions({
    collection: gridResult,
    reducer   : ee.Reducer.sum().setOutputs(['retama_' + year])
                  .combine(ee.Reducer.count().setOutputs(['total_' + year]), '', true),
    scale     : 10,
    tileScale : 4
  });
});

// =============================================================================
// 3. COMPUTE PROPORTIONS
// =============================================================================
exportYears.forEach(function(year) {
  gridResult = gridResult.map(function(cell) {
    var ret = ee.Number(cell.get('retama_' + year));
    var tot = ee.Number(cell.get('total_' + year));
    return cell.set('proportion_' + year, ret.divide(tot));
  });
});

// =============================================================================
// 4. VISUALISATION
// =============================================================================
var roi = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');
Map.centerObject(roi, 9);

// Show 2023 proportions as a choropleth
var vis2023 = gridResult.map(function(cell) {
  return cell.set('vis_proportion', cell.get('proportion_2023'));
});
Map.addLayer(
  vis2023,
  {color: 'red'},
  'Grid — invasion proportion 2023 (see Console)', false
);

print('Sample cell (first):', gridResult.first());

// =============================================================================
// 5. EXPORT
// =============================================================================
var exportName = '13_gridInvasion_' + CELL_SIZE + 'm' + RUN_SUFFIX;

Export.table.toAsset({
  collection : gridResult,
  description: exportName + '_asset',
  assetId    : ASSET_PREFIX + exportName
});

if (EXPORT_TO_DRIVE) {
  Export.table.toDrive({
    collection    : gridResult,
    description   : exportName + '_drive',
    folder        : DRIVE_FOLDER,
    fileNamePrefix: exportName,
    fileFormat    : 'SHP'
  });
}
