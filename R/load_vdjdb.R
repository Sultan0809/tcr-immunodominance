## R/load_vdjdb.R ------------------------------------------------------------
## Acquire and parse VDJdb, independent of which export the user supplied.
##
## VDJdb ships in several shapes:
##   * vdjdb.slim.txt        - built "slim" table
##   * vdjdb.txt             - built full table (has `meta`, `method`)
##   * SearchTable-*.tsv     - web browser export (different capitalisation)
## All three are handled. Columns are normalised to one canonical snake_case
## schema so nothing downstream has to care about the source.
## -----------------------------------------------------------------------------

## Canonical schema <- accepted aliases (post `normalise_names()`)
.VDJDB_ALIASES <- list(
  complex_id      = c("complex_id"),
  gene            = c("gene"),
  cdr3            = c("cdr3"),
  v_segm          = c("v_segm", "v", "v_segment"),
  j_segm          = c("j_segm", "j", "j_segment"),
  species         = c("species"),
  mhc_a           = c("mhc_a"),
  mhc_b           = c("mhc_b"),
  mhc_class       = c("mhc_class"),
  epitope         = c("antigen_epitope", "epitope"),
  epitope_gene    = c("antigen_gene", "epitope_gene"),
  epitope_species = c("antigen_species", "epitope_species"),
  reference_id    = c("reference_id", "reference"),
  method          = c("method", "web_method"),
  meta            = c("meta"),
  vdjdb_score     = c("vdjdb_score", "score")
)

#' Attempt to download VDJdb into data/raw
#'
#' Tries, in order: the GitHub release archive, then the repository master
#' archive. Both can fail (rate limits, network policy, or a release that ships
#' source requiring a build step). Failure is **not** an error — the function
#' returns FALSE and the caller falls back to a manual download, which is
#' documented in the README.
#'
#' @param dest_dir Directory to download into.
#' @param force Re-download even if a VDJdb file is already present.
#' @return Logical: is a usable VDJdb file now present in `dest_dir`?
download_vdjdb <- function(dest_dir = dir_raw(), force = FALSE) {
  ensure_dir(dest_dir)

  if (!force && length(find_vdjdb_file(dest_dir, quiet = TRUE))) {
    log_msg("VDJdb file already present in %s - skipping download", dest_dir)
    return(TRUE)
  }

  urls <- c(
    release = "https://github.com/antigenomics/vdjdb-db/releases/latest/download/vdjdb.zip",
    master  = "https://github.com/antigenomics/vdjdb-db/archive/refs/heads/master.zip"
  )

  zip_path <- file.path(dest_dir, "vdjdb_download.zip")
  for (nm in names(urls)) {
    log_msg("trying VDJdb source: %s", nm)
    if (try_download(urls[[nm]], zip_path)) {
      ok <- tryCatch({
        utils::unzip(zip_path, exdir = file.path(dest_dir, "vdjdb_unpacked"))
        TRUE
      }, error = function(e) {
        log_msg("unzip failed: %s", conditionMessage(e), level = "warn")
        FALSE
      })
      unlink(zip_path)
      if (ok && length(find_vdjdb_file(dest_dir, quiet = TRUE))) {
        log_msg("VDJdb acquired from %s", nm, level = "ok")
        return(TRUE)
      }
      log_msg(paste0(
        "%s downloaded, but it contains no assembled VDJdb table. The source ",
        "repository stores per-study chunks that require the VDJdb assembler ",
        "to build; the web export below is faster."), nm, level = "warn")
    }
  }

  log_msg(paste(
    "Automatic VDJdb download failed. Export it manually instead:",
    "  1. https://vdjdb.cdr3.net  ->  CDR3 tab",
    "  2. tick all boxes under 'General'",
    "  3. Export -> TSV",
    sprintf("  4. save SearchTable-*.tsv into %s", dest_dir),
    sep = "\n"
  ), level = "warn")
  FALSE
}

#' Locate a VDJdb table on disk
#'
#' Searches recursively and prefers the richest available format, because
#' `vdjdb.txt` carries the `meta` column (donor IDs) that the slim table drops.
#'
#' @param dir Directory to search.
#' @param quiet Suppress logging.
#' @return Character scalar path, or `character(0)` if none found.
#' Does this file actually look like a VDJdb table?
#'
#' Checks the header rather than the filename. The VDJdb repository contains
#' plenty of unrelated `.tsv` files (website assets, metadata) whose names would
#' otherwise match.
#'
#' @param path File to inspect.
#' @return Logical scalar.
is_vdjdb_table <- function(path) {
  if (!file.exists(path) || isTRUE(file.info(path)$isdir)) return(FALSE)
  hdr <- tryCatch(
    suppressWarnings(readLines(path, n = 1L, warn = FALSE)),
    error = function(e) character(0)
  )
  if (!length(hdr) || !nzchar(hdr[1])) return(FALSE)

  cols <- normalise_names(strsplit(hdr[1], "\t", fixed = TRUE)[[1]])
  has_cdr3    <- any(c("cdr3", "cdr3_beta", "cdr3_alpha") %in% cols)
  has_epitope <- any(c("antigen_epitope", "epitope") %in% cols)
  has_cdr3 && has_epitope
}

find_vdjdb_file <- function(dir = dir_raw(), quiet = FALSE) {
  if (!dir.exists(dir)) return(character(0))

  candidates <- list.files(
    dir, pattern = "\\.(txt|tsv)$",
    recursive = TRUE, full.names = TRUE, ignore.case = TRUE
  )
  ## Website assets and documentation live alongside the data in the source repo.
  candidates <- candidates[!grepl("/(images|docs|\\.github)/", candidates)]
  candidates <- candidates[file.info(candidates)$size > 1000]
  if (!length(candidates)) return(character(0))

  ## Validate by header, not by name.
  valid <- candidates[vapply(candidates, is_vdjdb_table, logical(1))]
  if (!length(valid)) {
    if (!quiet) {
      log_msg("found %d .tsv/.txt files under %s, but none has VDJdb columns",
              length(candidates), dir, level = "warn")
    }
    return(character(0))
  }

  ## Preference: built full table > web export > slim > largest remaining
  ## (chunk files are small; a concatenated database is large).
  rank_of <- function(p) {
    b <- tolower(basename(p))
    if (identical(b, "vdjdb.txt"))      return(1L)
    if (grepl("^searchtable", b))        return(2L)
    if (grepl("full", b))                return(3L)
    if (grepl("slim", b))                return(4L)
    5L
  }
  ord <- order(vapply(valid, rank_of, integer(1)),
               -file.info(valid)$size)
  chosen <- valid[ord][1]
  if (!quiet) {
    log_msg("using VDJdb file: %s (%.1f MB, %d valid candidates)",
            chosen, file.info(chosen)$size / 1e6, length(valid))
  }
  chosen
}

#' Extract a field from VDJdb's JSON `meta` column
#'
#' Uses regex rather than a per-row JSON parse: the `meta` column has ~100k rows
#' and is frequently malformed or empty, which makes `fromJSON()` both slow and
#' fragile here.
#'
#' @param meta Character vector of JSON strings.
#' @param field Field name as it appears in the JSON, e.g. "subject.id".
#' @return Character vector, `NA` where absent or empty.
extract_meta_field <- function(meta, field) {
  if (is.null(meta)) return(character(0))
  meta <- as.character(meta)
  out <- rep(NA_character_, length(meta))
  if (!length(meta)) return(out)

  field_rx <- gsub("\\.", "\\\\.", field)
  pattern  <- sprintf('"%s"\\s*:\\s*"([^"]*)"', field_rx)

  ## `regmatches()` drops non-matching elements, which would misalign the
  ## output. Masking with grepl() and substituting in place keeps alignment.
  has <- !is.na(meta) & grepl(pattern, meta, perl = TRUE)
  if (any(has)) {
    out[has] <- sub(paste0("^.*?", pattern, ".*$"), "\\1", meta[has], perl = TRUE)
  }
  out[!is.na(out) & !nzchar(out)] <- NA_character_
  out
}

#' Load and normalise VDJdb
#'
#' @param path Path to a VDJdb table. Defaults to whatever `find_vdjdb_file()`
#'   turns up in `data/raw`.
#' @param species Keep only this species (default "HomoSapiens"). `NULL` keeps all.
#' @param chain Keep only this chain, "TRA" or "TRB". `NULL` keeps both.
#' @param min_score Minimum `vdjdb_score` to retain. VDJdb scores run 0-3, where
#'   0 means "critical information missing". Default 1 drops the lowest tier.
#' @return A `data.table` with the canonical schema plus parsed meta columns.
load_vdjdb <- function(path = NULL,
                       species = "HomoSapiens",
                       chain = "TRB",
                       min_score = 1L) {

  if (is.null(path)) path <- find_vdjdb_file()
  if (!length(path) || !file.exists(path)) {
    stop("No VDJdb file found. Run download_vdjdb(), or export one manually ",
         "(see README). Expected somewhere under data/raw/.", call. = FALSE)
  }

  raw <- data.table::fread(path, sep = "\t", quote = "", na.strings = c("", "NA"))
  data.table::setnames(raw, normalise_names(names(raw)))
  log_msg("read %d raw rows, %d columns", nrow(raw), ncol(raw))

  ## --- map aliases onto the canonical schema -------------------------------
  cols <- lapply(names(.VDJDB_ALIASES), function(canon) {
    hit <- intersect(.VDJDB_ALIASES[[canon]], names(raw))
    if (length(hit)) as.character(raw[[hit[1]]]) else rep(NA_character_, nrow(raw))
  })
  names(cols) <- names(.VDJDB_ALIASES)
  out <- data.table::as.data.table(cols)

  ## A required column that is absent, or present but entirely empty, means the
  ## file is not a usable VDJdb export.
  core <- c("cdr3", "epitope", "gene")
  empty_core <- core[vapply(core, function(cn) all(is.na(out[[cn]])), logical(1))]
  if (length(empty_core)) {
    stop("VDJdb file has no usable values in required column(s): ",
         paste(empty_core, collapse = ", "),
         ". The file may be truncated or not a VDJdb export.", call. = FALSE)
  }

  ## --- types ----------------------------------------------------------------
  out[, vdjdb_score := suppressWarnings(as.integer(vdjdb_score))]
  out[, complex_id  := as.character(complex_id)]

  ## --- parse donor / study metadata ----------------------------------------
  ## `subject_id` is the basis for every publicness metric in step 2. It is only
  ## unique *within* a study, so the usable donor key is study + subject.
  out[, subject_id := extract_meta_field(meta, "subject.id")]
  out[, study_id   := extract_meta_field(meta, "study.id")]
  out[, tissue     := extract_meta_field(meta, "tissue")]
  out[, donor_mhc  := extract_meta_field(meta, "donor.MHC")]

  out[, has_subject_id := !is.na(subject_id)]
  out[, donor_key := data.table::fifelse(
    has_subject_id,
    paste(data.table::fcoalesce(study_id, reference_id, "NA"), subject_id, sep = "|"),
    NA_character_
  )]

  ## --- filters --------------------------------------------------------------
  ## NB: `species` is both an argument and a column name. data.table's `..`
  ## prefix only works in `j`, not `i`, so the argument is copied to a
  ## differently-named local to avoid the column masking it.
  n0 <- nrow(out)
  keep_species <- species
  keep_chain   <- chain
  if (!is.null(keep_species)) out <- out[species == keep_species]
  if (!is.null(keep_chain))   out <- out[gene == keep_chain]
  if (!is.null(min_score))    out <- out[!is.na(vdjdb_score) & vdjdb_score >= min_score]

  ## Drop rows unusable for any downstream analysis.
  out <- out[!is.na(cdr3) & nzchar(cdr3) & !is.na(epitope) & nzchar(epitope)]

  ## Guard against non-amino-acid CDR3s (a small number of malformed records).
  out <- out[grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", cdr3)]

  log_msg("after filtering: %d rows (%.1f%% of raw)",
          nrow(out), 100 * nrow(out) / max(n0, 1L))
  log_msg("  epitopes: %d | donors with IDs: %d | studies: %d",
          data.table::uniqueN(out$epitope),
          data.table::uniqueN(out$donor_key, na.rm = TRUE),
          data.table::uniqueN(out$reference_id))

  out[]
}
