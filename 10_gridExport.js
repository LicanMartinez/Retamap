// =============================================================================
// =============================================================================
// 10: Generate and export a regular grid over the study area
// =============================================================================
// -----------------------------------------------------------------------------
// Produces a vector grid (FeatureCollection) covering the bounding box of
// '3_study_area_retama', computed in UTM 19S (EPSG:32719).
//
// Output:
//   Asset  → ASSET_PREFIX + 'grid_' + CELL_SIZE + 'm'
//   Drive  → DRIVE_FOLDER + '/grid_' + CELL_SIZE + 'm'  (GeoJSON)
//
// Usage: set CELL_SIZE to 1000, 3000, or 7000 (metres) and run.
// The cell size is used as an export suffix so multiple runs coexist.
// =============================================================================


// =============================================================================
// 0. CONFIGURATION
// =============================================================================
var ASSET_PREFIX = 'projects/ee-licanemartinez/assets/Retamap/';
var DRIVE_FOLDER = 'Retamap/Retamap_GEE_Exports';

// ── Edit this value before running ──────────────────────────────────────────
var CELL_SIZE = 3000;   // metres — try 1000, 3000, 7000
// ─────────────────────────────────────────────────────────────────────────────

var UTM19S = ee.Projection('EPSG:32719');

var roi = ee.FeatureCollection(ASSET_PREFIX + '3_study_area_retama');

// =============================================================================
// 1. BOUNDING BOX IN UTM 19S
// =============================================================================
// .bounds() in the target projection returns the axis-aligned bbox.
// errorMargin of 1 m is sufficient for a bounding-box operation.
var bboxUTM = roi.geometry().bounds(1, UTM19S);

// =============================================================================
// 2. BUILD GRID
// =============================================================================
// coveringGrid() tiles the geometry with cells of size CELL_SIZE in the
// given projection.  Each Feature carries a 'grid_id' property (col_row).
var gridRaw = bboxUTM.coveringGrid(UTM19S, CELL_SIZE);

// Add an integer id and the cell size as metadata to every cell
var grid = gridRaw.map(function(cell) {
  return cell.set('cell_size_m', CELL_SIZE);
}).map(function(cell, idx) {
  // sequential id via iterate alternative: use system:index (already present)
  return cell.set('grid_id', cell.get('system:index'));
});

// =============================================================================
// 3. VISUALISATION
// =============================================================================
Map.centerObject(roi, 9);
Map.addLayer(roi,  {color: 'grey'},  'Study area', false);
Map.addLayer(grid, {color: '1a9641', fillColor: '00000000'}, 'Grid ' + CELL_SIZE + ' m');

print('Grid cell count:', grid.size());
print('Projection:', UTM19S.wkt());
print('BBox (UTM 19S):', bboxUTM);

// =============================================================================
// 4. EXPORT
// =============================================================================
var exportName = 'grid_' + CELL_SIZE + 'm';

// ── 4a. Asset (FeatureCollection) ────────────────────────────────────────────
Export.table.toAsset({
  collection: grid,
  description: exportName + '_asset',
  assetId:     ASSET_PREFIX + exportName
});

// ── 4b. Drive (GeoJSON) ──────────────────────────────────────────────────────
Export.table.toDrive({
  collection:  grid,
  description: exportName + '_drive',
  folder:      DRIVE_FOLDER,
  fileNamePrefix: exportName,
  fileFormat:  'GeoJSON'
});
