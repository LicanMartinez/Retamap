# =============================================================================
# 09v_valRandomize.R
# -----------------------------------------------------------------------------
# Step 2 of the independent-validation workflow. Consumes the raw point CSV
# exported by 09v_valGenerator.js (Google Drive) and produces:
#
#   1. 09v_master_key.csv        — pxid + stratum + map labels (m2017,m2025) +
#                                  lon/lat + assigned validators. LOCAL ONLY,
#                                  never uploaded/published (it carries labels).
#   2. fc/09v_valPoints_fc.shp   — BLIND FeatureCollection: pixel-contour
#      (+ .geojson)                polygons with ONLY `pxid`. Upload as GEE asset
#                                  `09v_valPoints_fc`. Centroid is derived in GEE
#                                  with .centroid().
#   3. sheets/09v_sheet_<val>.csv — one CSV per validator → one tab of the
#                                  Google Sheet. Columns to fill:
#                                  es_borde, class_2017, class_2017_confidence,
#                                  class_2025, class_2025_confidence.
#   4. kml/<validator>/<pxid>.kml — per-point KML (centroid + 10 m contour) for
#                                  Google Earth Pro historical imagery.
#
# Design: 3 strata, N drawn by Olofsson (set in the generator). 3 validators;
# each point is assigned to 2 of the 3 (round-robin pairs AB/BC/CA), balanced.
# pxid is a SHUFFLED sequential id so strata are not contiguous in the sheets.
# =============================================================================

suppressPackageStartupMessages(library(sf))

# =============================================================================
# 0. CONFIGURATION  (edit these)
# =============================================================================
RUN_SUFFIX <- "_trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k"

PROJECT_DIR <- "D:/Lican/uni/Investigacion/Colaboraciones_y_ayudas/Retamap"
WORK_DIR    <- file.path(PROJECT_DIR, "gee_scripts", "validation_independent")

# Place the Drive-downloaded raw CSV here (default name from the generator):
INPUT_CSV <- file.path(WORK_DIR, paste0("09v_valPoints_raw", RUN_SUFFIX, ".csv"))

VALIDATORS <- c("lican", "sofi", "val3")   # 3 names → 3 sheets / 3 kml folders
SEED       <- 42
PIXEL_SIZE <- 10                            # m, Sentinel-2 pixel side
UTM_EPSG   <- 32719                         # UTM 19S

# =============================================================================
# 1. READ RAW POINTS
# =============================================================================
stopifnot(file.exists(INPUT_CSV))
raw <- read.csv(INPUT_CSV, stringsAsFactors = FALSE)
# Expected columns: stratum, m2017, m2025, lon, lat (+ possibly system:index/.geo)
need <- c("stratum", "lon", "lat")
stopifnot(all(need %in% names(raw)))
for (cc in c("m2017", "m2025")) if (!cc %in% names(raw)) raw[[cc]] <- NA
raw <- raw[, c("stratum", "m2017", "m2025", "lon", "lat")]
cat(sprintf("Read %d points (S1=%d, S2=%d, S3=%d)\n",
            nrow(raw), sum(raw$stratum == 1), sum(raw$stratum == 2), sum(raw$stratum == 3)))

# =============================================================================
# 2. SHUFFLE → pxid (strata not contiguous)
# =============================================================================
set.seed(SEED)
ord <- sample(seq_len(nrow(raw)))
raw <- raw[ord, ]
width <- nchar(as.character(nrow(raw)))
raw$pxid <- sprintf(paste0("px%0", width, "d"), seq_len(nrow(raw)))
rownames(raw) <- NULL

# =============================================================================
# 3. ASSIGN VALIDATORS (each point to 2 of 3, round-robin pairs)
# =============================================================================
pairs <- list(c(1, 2), c(2, 3), c(1, 3))             # AB, BC, CA
pidx  <- (seq_len(nrow(raw)) - 1) %% 3 + 1
raw$validators <- vapply(pidx, function(k)
  paste(VALIDATORS[pairs[[k]]], collapse = ";"), character(1))
cat("Validator load (points each):\n")
print(sapply(VALIDATORS, function(v) sum(grepl(v, raw$validators))))

# =============================================================================
# 4. GEOMETRY: centroid (4326) + 10 m pixel contour (built in UTM)
# =============================================================================
pts_utm <- st_transform(st_as_sf(raw, coords = c("lon", "lat"), crs = 4326), UTM_EPSG)
co      <- st_coordinates(pts_utm)
half    <- PIXEL_SIZE / 2

polys <- lapply(seq_len(nrow(co)), function(i) {
  x <- co[i, 1]; y <- co[i, 2]
  st_polygon(list(matrix(c(
    x - half, y - half,  x + half, y - half,
    x + half, y + half,  x - half, y + half,
    x - half, y - half), ncol = 2, byrow = TRUE)))
})
contour_utm  <- st_sf(pxid = raw$pxid, geometry = st_sfc(polys, crs = UTM_EPSG))
contour_4326 <- st_transform(contour_utm, 4326)

# =============================================================================
# 5. WRITE OUTPUTS
# =============================================================================
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(WORK_DIR, "fc"),     showWarnings = FALSE)
dir.create(file.path(WORK_DIR, "sheets"), showWarnings = FALSE)
dir.create(file.path(WORK_DIR, "kml"),    showWarnings = FALSE)

# 5a. Master key (LOCAL ONLY — carries labels) -------------------------------
write.csv(raw[, c("pxid", "stratum", "m2017", "m2025", "lon", "lat", "validators")],
          file.path(WORK_DIR, "09v_master_key.csv"), row.names = FALSE, na = "")

# 5b. Blind FeatureCollection (pxid only) → asset upload -----------------------
fc_shp <- file.path(WORK_DIR, "fc", "09v_valPoints_fc.shp")
fc_gj  <- file.path(WORK_DIR, "fc", "09v_valPoints_fc.geojson")
st_write(contour_4326, fc_shp, delete_dsn = TRUE, quiet = TRUE)
st_write(contour_4326, fc_gj,  delete_dsn = TRUE, quiet = TRUE)

# 5c. Per-validator sheets (blind; columns to fill) ---------------------------
empty_cols <- c("es_borde", "class_2017", "class_2017_confidence",
                "class_2025", "class_2025_confidence")
for (v in VALIDATORS) {
  sub <- raw[grepl(v, raw$validators), "pxid", drop = FALSE]
  for (cc in empty_cols) sub[[cc]] <- NA
  write.csv(sub, file.path(WORK_DIR, "sheets", paste0("09v_sheet_", v, ".csv")),
            row.names = FALSE, na = "")
}

# 5d. Per-point KML (centroid + contour) in each assigned validator's folder ---
ring_all <- st_coordinates(contour_4326)   # cols X,Y,...,L2 (feature index)
fid_col  <- ncol(ring_all)                 # last col = feature id (L2)
make_kml <- function(pxid, cx, cy, ring) {
  cs <- paste(apply(ring, 1, function(r) paste0(r[1], ",", r[2], ",0")), collapse = " ")
  paste0(
    '<?xml version="1.0" encoding="UTF-8"?>\n',
    '<kml xmlns="http://www.opengis.net/kml/2.2"><Document><name>', pxid, '</name>\n',
    '<Placemark><name>', pxid, ' centroid</name>',
    '<Point><coordinates>', cx, ',', cy, ',0</coordinates></Point></Placemark>\n',
    '<Placemark><name>', pxid, ' pixel</name>',
    '<Style><LineStyle><color>ff0000ff</color><width>2</width></LineStyle>',
    '<PolyStyle><fill>0</fill></PolyStyle></Style>',
    '<Polygon><outerBoundaryIs><LinearRing><coordinates>', cs,
    '</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>\n',
    '</Document></kml>')
}
for (v in VALIDATORS) dir.create(file.path(WORK_DIR, "kml", v), showWarnings = FALSE)
for (i in seq_len(nrow(raw))) {
  ring <- ring_all[ring_all[, fid_col] == i, c("X", "Y"), drop = FALSE]
  kml  <- make_kml(raw$pxid[i], raw$lon[i], raw$lat[i], ring)
  for (v in VALIDATORS[pairs[[pidx[i]]]]) {
    writeLines(kml, file.path(WORK_DIR, "kml", v, paste0(raw$pxid[i], ".kml")))
  }
}

cat(sprintf("\nDone. Outputs in: %s\n", WORK_DIR))
cat("  - 09v_master_key.csv  (LOCAL ONLY, has labels)\n")
cat("  - fc/09v_valPoints_fc.{shp,geojson}  → upload as GEE asset 09v_valPoints_fc\n")
cat(sprintf("  - sheets/09v_sheet_<%s>.csv  → import as tabs of one Google Sheet\n",
            paste(VALIDATORS, collapse = "/")))
cat("  - kml/<validator>/<pxid>.kml  → open in Google Earth Pro\n")
