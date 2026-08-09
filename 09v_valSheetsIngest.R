# =============================================================================
# 09v_valSheetsIngest.R
# -----------------------------------------------------------------------------
# Turns the validator Google Sheets, pasted verbatim into text files, into the
# canonical `<validator>_filled.csv` that 09v_valAnalysis_compute.R consumes.
#
#   validation/sheets_filled/<val>_pasted.txt   (in:  tab-separated paste)
#   validation/sheets_filled/<val>_filled.csv   (out: validator,pxid,
#                                                class_2017,conf_2017,
#                                                class_2025,conf_2025)
#
# Why this exists: there is no automated download path. The Drive MCP connector
# is not present in every session, and even when it is, its CSV export returns
# the wrong tab of a multi-tab Sheet and its markdown reader truncates large
# sheets silently (Sofi's ~540 rows). Pasting the data tab into a text file and
# running this script sidesteps all of that AND reconciles the result against
# the master key, so a bad paste fails loudly instead of quietly shrinking n.
#
# Handles the real quirks of the current sheets:
#   - CRLF line endings and a possible UTF-8 BOM from the browser paste
#   - `to_sofi` column in Lican's and Jaime's sheets: rows marked "no validar"
#     were reassigned to Sofi on 2026-07-27 and must NOT be counted as theirs
#   - `observaciones` free-text column present in some sheets, absent in others
#   - duplicated pxid rows (a paste artefact): collapsed when the duplicates are
#     identical, hard error when they disagree
#   - rows not yet labelled: kept out of the output (the engine drops them too)
#
# Run:  Rscript gee_scripts/09v_valSheetsIngest.R
# =============================================================================

PROJECT_DIR <- "D:/Lican/uni/Investigacion/Colaboraciones_y_ayudas/Retamap"
WORK_DIR    <- file.path(PROJECT_DIR, "gee_scripts", "validation")
SHEETS_DIR  <- file.path(WORK_DIR, "sheets_filled")
MASTER_KEY  <- file.path(WORK_DIR, "09v_master_key.csv")
VALIDATORS  <- c("lican", "sofi", "jaime")
SKIP_FLAG   <- "no validar"     # value of `to_sofi` marking a reassigned row

mk <- read.csv(MASTER_KEY, stringsAsFactors = FALSE, colClasses = "character")
mk_val <- setNames(strsplit(mk$validators, ";"), mk$pxid)

trim <- function(x) gsub("^\\s+|\\s+$", "", x)

read_pasted <- function(path) {
  txt <- readLines(path, warn = FALSE, encoding = "UTF-8")
  txt <- sub("^\ufeff", "", txt)          # strip BOM if the paste carried one
  txt <- sub("\r$", "", txt)              # CRLF -> LF
  txt <- txt[trim(txt) != ""]
  parts <- strsplit(txt, "\t", fixed = TRUE)
  hdr   <- trim(parts[[1]])
  # pad short rows: trailing empty cells get dropped by the paste
  body  <- lapply(parts[-1], function(p) {
    p <- trim(p); length(p) <- length(hdr); ifelse(is.na(p), "", p) })
  df <- as.data.frame(do.call(rbind, body), stringsAsFactors = FALSE)
  names(df) <- hdr
  df
}

summary_rows <- list()
for (v in VALIDATORS) {
  src <- file.path(SHEETS_DIR, paste0(v, "_pasted.txt"))
  if (!file.exists(src)) { message("no pasted file for ", v, " — skipping"); next }
  df <- read_pasted(src)
  need <- c("pxid", "class_2017", "class_2017_confidence",
            "class_2025", "class_2025_confidence")
  miss <- setdiff(need, names(df))
  if (length(miss)) stop(v, ": pasted sheet is missing column(s): ",
                         paste(miss, collapse = ", "))
  n_raw <- nrow(df)

  # 1. drop rows reassigned to another validator
  n_reassigned <- 0L
  if ("to_sofi" %in% names(df)) {
    drop <- df$to_sofi == SKIP_FLAG
    n_reassigned <- sum(drop)
    df <- df[!drop, , drop = FALSE]
  }

  # 2. collapse duplicated pxids, but only when they carry the same labels
  key <- df[, need]
  dups <- unique(df$pxid[duplicated(df$pxid)])
  for (p in dups) {
    rows <- unique(key[key$pxid == p, ])
    if (nrow(rows) > 1)
      stop(v, ": pxid ", p, " appears more than once with DIFFERENT labels — ",
           "resolve it in the sheet before re-running")
  }
  n_dups <- sum(duplicated(df$pxid))
  df <- df[!duplicated(df$pxid), , drop = FALSE]

  # 3. reconcile against the master key (authoritative assignment)
  assigned <- names(mk_val)[vapply(mk_val, function(a) v %in% a, logical(1))]
  extra   <- setdiff(df$pxid, assigned)
  missing <- setdiff(assigned, df$pxid)
  if (length(extra))
    stop(v, ": ", length(extra), " pxid(s) in the sheet are NOT assigned to ",
         v, " in the master key: ", paste(head(extra, 10), collapse = ", "))
  if (length(missing))
    warning(v, ": ", length(missing), " assigned pxid(s) absent from the sheet: ",
            paste(head(missing, 10), collapse = ", "))

  # 4. keep only rows with at least one year labelled
  has <- trim(df$class_2017) != "" | trim(df$class_2025) != ""
  out <- data.frame(
    validator   = v,
    pxid        = df$pxid[has],
    class_2017  = trim(df$class_2017[has]),
    conf_2017   = trim(df$class_2017_confidence[has]),
    class_2025  = trim(df$class_2025[has]),
    conf_2025   = trim(df$class_2025_confidence[has]),
    stringsAsFactors = FALSE)

  bad <- setdiff(unique(c(out$class_2017, out$class_2025)),
                 c("Retama", "No_retama", ""))
  if (length(bad)) stop(v, ": unexpected class value(s): ",
                        paste(bad, collapse = ", "))

  dst <- file.path(SHEETS_DIR, paste0(v, "_filled.csv"))
  write.csv(out, dst, row.names = FALSE)
  summary_rows[[v]] <- data.frame(
    validator = v, rows_pasted = n_raw, reassigned_dropped = n_reassigned,
    duplicates_collapsed = n_dups, assigned_in_master_key = length(assigned),
    still_unlabelled = length(assigned) - nrow(out),
    written = nrow(out),
    pct_complete = round(100 * nrow(out) / length(assigned), 1),
    row.names = NULL)
}

res <- do.call(rbind, summary_rows)
cat("\nIngest summary\n"); print(res, row.names = FALSE)
cat("\nWrote <validator>_filled.csv to:", SHEETS_DIR, "\n")
