// =============================================================================
// =============================================================================
// 12: Generate and export a regular grid over the study area
// =============================================================================
// -----------------------------------------------------------------------------
// Produces a vector grid (FeatureCollection) covering the bounding box of
// '3_study_area_retama', computed in UTM 19S (EPSG:32719).
//
// CLIP_TO_ROI = false  →  full rectangular cells over the bbox
// CLIP_TO_ROI = true   →  cells clipped to the exact ROI polygon shape;
//                          boundary cells are trimmed, tiny slivers are dropped
//
// Output:
//   Asset  → ASSET_PREFIX + '12_grid_' + CELL_SIZE + 'm'         (bbox)
//          or ASSET_PREFIX + '12_grid_' + CELL_SIZE + 'm_clipped' (clipped)
//   Drive  → same naming, GeoJSON (if EXPORT_TO_DRIVE = true)
//
// Usage: set CELL_SIZE to 1000, 3000, or 7000 (metres) and run.
// =============================================================================


// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX    = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER    = 'Retamap/Retamap_GEE_Exports';
var EXPORT_TO_DRIVE = false;  // toggle: true → also export to Google Drive

// ── Edit these values before running ─────────────────────────────────────────
var CELL_SIZE   = 3000;   // metres — try 1000, 3000, 7000
var CLIP_TO_ROI = true;   // true → clip cells to the actual ROI polygon shape
// ─────────────────────────────────────────────────────────────────────────────

var UTM19S = ee.Projection('EPSG:32719');

var roi     = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');
var roiGeom = roi.geometry();   // dissolved geometry (may be multi-polygon)

// =============================================================================
// 1. BOUNDING BOX IN UTM 19S
// =============================================================================
var bboxUTM = roiGeom.bounds(1, UTM19S);

// =============================================================================
// 2. BUILD GRID (rectangular, over bbox)
// =============================================================================
var gridRaw = bboxUTM.coveringGrid(UTM19S, CELL_SIZE);

// Add metadata to every cell
var gridBbox = gridRaw.map(function(cell) {
  return cell
    .set('cell_size_m', CELL_SIZE)
    .set('grid_id', cell.get('system:index'));
});

// =============================================================================
// 3. OPTIONAL CLIP TO ROI SHAPE
// =============================================================================
var grid;

if (CLIP_TO_ROI) {
  var MIN_AREA_M2 = CELL_SIZE * CELL_SIZE * 0.01;

  grid = gridBbox
    .filterBounds(roiGeom)
    .map(function(cell) {
      var clippedGeom = cell.geometry().intersection(roiGeom, 1);
      return ee.Feature(clippedGeom)
        .copyProperties(cell)
        .set('area_m2', clippedGeom.area(1));
    })
    .filter(ee.Filter.gt('area_m2', MIN_AREA_M2));

} else {
  grid = gridBbox;
}

// =============================================================================
// 4. VISUALISATION
// =============================================================================
Map.centerObject(roi, 9);
Map.addLayer(roi,      {color: '555555'},                        'Study area ROI', true);
Map.addLayer(gridBbox, {color: 'aaaaaa', fillColor: '00000000'}, 'Grid bbox ' + CELL_SIZE + ' m', !CLIP_TO_ROI);
if (CLIP_TO_ROI) {
  Map.addLayer(grid, {color: '1a9641', fillColor: '1a964122'},   'Grid clipped ' + CELL_SIZE + ' m', true);
}

print('Bbox cell count:',    gridBbox.size());
if (CLIP_TO_ROI) {
  print('Clipped cell count:', grid.size());
}

// =============================================================================
// 5. EXPORT (asset name prefix: 12_)
// =============================================================================
var clipSuffix = CLIP_TO_ROI ? '_clipped' : '';
var exportName = '12_grid_' + CELL_SIZE + 'm' + clipSuffix;

// ── 5a. Asset (FeatureCollection) ────────────────────────────────────────────
Export.table.toAsset({
  collection : grid,
  description: exportName + '_asset',
  assetId    : ASSET_PREFIX + exportName
});

// ── 5b. Drive (GeoJSON) ──────────────────────────────────────────────────────
if (EXPORT_TO_DRIVE) {
  Export.table.toDrive({
    collection    : grid,
    description   : exportName + '_drive',
    folder        : DRIVE_FOLDER,
    fileNamePrefix: exportName,
    fileFormat    : 'GeoJSON'
  });
}
