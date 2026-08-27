# =============================================================================
# 09v_valKmlRefresh.R
# -----------------------------------------------------------------------------
# Regenerates ONLY the per-point KMLs (exact grid-A / 4326 prediction pixel)
# from a fresh raw CSV that adds lonA/latA/dLon/dLat, WITHOUT re-shuffling or
# re-assigning pxids. The existing pxid <-> validator assignment
# (09v_master_key.csv, already handed out to lican/sofi/jaime) must not change
# — see 09v_valRandomize.R header. This script only reads that assignment and
# writes fresh KML geometry.
#
# Matches rows between the master key and the new raw CSV by (stratum, lon,
# lat) — those columns are untouched by the generator re-run, only
# lonA/latA/dLon/dLat are new. Output goes into kml/<validator>_new/ so the
# original kml/<validator>/ folders are left untouched for comparison.
# =============================================================================

# Folder that CONTAINS gee_scripts/. Defaults to the directory the script is
# invoked from (the project root); override with the RETAMAP_DIR env variable.
PROJECT_DIR <- Sys.getenv("RETAMAP_DIR", unset = ".")
WORK_DIR    <- file.path(PROJECT_DIR, "gee_scripts", "validation_data")
RUN_SUFFIX  <- "_trimmedRF1comp_s1n10k_qs1n10k_s0n15k_qs0n15k"

MASTER_KEY_CSV <- file.path(WORK_DIR, "09v_master_key.csv")
NEW_RAW_CSV    <- file.path(WORK_DIR, paste0("09v_valPoints_raw", RUN_SUFFIX, "_new.csv"))
VALIDATORS     <- c("lican", "sofi", "jaime")

# =============================================================================
# 1. LOAD + JOIN (no shuffle, no re-assignment)
# =============================================================================
stopifnot(file.exists(MASTER_KEY_CSV), file.exists(NEW_RAW_CSV))
mk  <- read.csv(MASTER_KEY_CSV, stringsAsFactors = FALSE)
new <- read.csv(NEW_RAW_CSV,    stringsAsFactors = FALSE)

need <- c("stratum", "lon", "lat", "lonA", "latA", "dLon", "dLat")
stopifnot(all(need %in% names(new)))   # re-run 09v_valGenerator.js if lonA.. missing

# Join key: round lon/lat to 8 decimals (~1 mm) to absorb float round-trip noise.
rnd <- function(x) round(x, 8)
mk$key  <- paste(mk$stratum,  rnd(mk$lon),  rnd(mk$lat))
new$key <- paste(new$stratum, rnd(new$lon), rnd(new$lat))

stopifnot(!anyDuplicated(mk$key), !anyDuplicated(new$key))
m <- match(mk$key, new$key)
if (anyNA(m)) {
  stop(sprintf("%d/%d master-key points did NOT match the new CSV by (stratum,lon,lat) — check RUN_SUFFIX / that the new export used the same point set.",
               sum(is.na(m)), nrow(mk)))
}
cat(sprintf("Matched all %d master-key points to the new CSV.\n", nrow(mk)))

merged <- cbind(mk[, c("pxid", "validators")], new[m, c("lonA", "latA", "dLon", "dLat")])

# =============================================================================
# 2. KML BUILDERS (identical to 09v_valRandomize.R section 5d)
# =============================================================================
make_kml <- function(pxid, cx, cy, ring) {
  cs <- paste(apply(ring, 1, function(r) paste0(r[1], ",", r[2], ",0")), collapse = " ")
  paste0(
    '<?xml version="1.0" encoding="UTF-8"?>\n',
    '<kml xmlns="http://www.opengis.net/kml/2.2"><Document><name>', pxid, '</name>\n',
    '<Placemark><name>', pxid, ' center</name>',
    '<Point><coordinates>', cx, ',', cy, ',0</coordinates></Point></Placemark>\n',
    '<Placemark><name>', pxid, ' pixel</name>',
    '<Style><LineStyle><color>ff0000ff</color><width>2</width></LineStyle>',
    '<PolyStyle><fill>0</fill></PolyStyle></Style>',
    '<Polygon><outerBoundaryIs><LinearRing><coordinates>', cs,
    '</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>\n',
    '</Document></kml>')
}

# Axis-aligned (lon/lat) rectangle of the grid-A pixel around (lonA, latA).
gridA_ring <- function(lonA, latA, dLon, dLat) {
  hx <- dLon / 2; hy <- dLat / 2
  matrix(c(lonA - hx, latA - hy,   lonA + hx, latA - hy,
           lonA + hx, latA + hy,   lonA - hx, latA + hy,
           lonA - hx, latA - hy), ncol = 2, byrow = TRUE)
}

# =============================================================================
# 3. WRITE  kml/<validator>_new/<pxid>.kml
# =============================================================================
for (v in VALIDATORS) {
  dir.create(file.path(WORK_DIR, "kml", paste0(v, "_new")),
             recursive = TRUE, showWarnings = FALSE)
}

for (i in seq_len(nrow(merged))) {
  ring <- gridA_ring(merged$lonA[i], merged$latA[i], merged$dLon[i], merged$dLat[i])
  kml  <- make_kml(merged$pxid[i], merged$lonA[i], merged$latA[i], ring)
  for (v in strsplit(merged$validators[i], ";")[[1]]) {
    writeLines(kml, file.path(WORK_DIR, "kml", paste0(v, "_new"), paste0(merged$pxid[i], ".kml")))
  }
}

cat(sprintf("\nDone. %d KMLs written into kml/<%s>_new/\n",
            nrow(merged), paste(VALIDATORS, collapse = ">_new/, <")))
