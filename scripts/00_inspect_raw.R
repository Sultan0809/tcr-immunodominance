## scripts/00_inspect_raw.R ---------------------------------------------------
## Diagnostic: what is actually in data/raw, and which files look like VDJdb?
## Run this when the loader cannot find a usable table.
##
## Run:  Rscript scripts/00_inspect_raw.R
## -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

source(file.path("R", "utils.R"))
source(file.path("R", "load_vdjdb.R"))

log_msg("=== 00_inspect_raw ===")

root <- dir_raw()
if (!dir.exists(root)) {
  stop("data/raw does not exist - run scripts/01_download_data.R first.",
       call. = FALSE)
}

files <- list.files(root, pattern = "\\.(txt|tsv)$",
                    recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

if (!length(files)) {
  message("No .txt/.tsv files under ", root)
} else {
  info <- data.table(
    path     = sub(paste0("^", root, "/?"), "", files),
    size_kb  = round(file.info(files)$size / 1024, 1),
    is_vdjdb = vapply(files, is_vdjdb_table, logical(1))
  )
  setorder(info, -is_vdjdb, -size_kb)

  message("\n--- 20 largest .tsv/.txt files under data/raw ---")
  print(head(info, 20))

  message("\nFiles with valid VDJdb headers: ", sum(info$is_vdjdb))

  if (any(info$is_vdjdb)) {
    best <- files[vapply(files, is_vdjdb_table, logical(1))][1]
    message("\n--- header of first valid file ---")
    message(basename(best))
    print(normalise_names(strsplit(readLines(best, n = 1L, warn = FALSE),
                                   "\t", fixed = TRUE)[[1]]))
  } else {
    message("\nNo assembled VDJdb table present.")
    message("Get one from the web browser instead:")
    message("  1. https://vdjdb.cdr3.net  ->  CDR3 tab")
    message("  2. tick every box under 'General'")
    message("  3. Export -> TSV")
    message("  4. save the SearchTable-*.tsv into ", root)
    message("  5. re-run: Rscript scripts/02_parse_vdjdb.R")
  }
}

## Top-level directory listing, for orientation.
message("\n--- top-level contents of data/raw ---")
print(list.files(root, recursive = FALSE))
