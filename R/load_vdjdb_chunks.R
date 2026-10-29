## R/load_vdjdb_chunks.R -----------------------------------------------------
## Read VDJdb's per-study submission chunks directly.
##
## The VDJdb source repository stores one TSV per publication under `chunks/`,
## plus verified non-binders under `chunks_negative/`. Assembling them with the
## official tool needs a Python toolchain; reading them directly does not, and
## the chunk format is in several ways richer than the assembled table:
##
##   * `meta_*` fields are real columns, not a JSON blob to regex
##   * alpha and beta chains are paired on one row
##   * `chunks_negative/` exposes verified non-binders
##
## What the chunks do NOT carry is `vdjdb.score`, which the assembler derives
## from the `method_*` fields. Confidence filtering therefore works off those
## fields directly - see `summarise_methods()` and `chunk_confidence()`.
## -----------------------------------------------------------------------------

## Chunk column -> canonical name. Anything unlisted is dropped.
.CHUNK_MAP <- c(
  species               = "species",
  mhc_a                 = "mhc_a",
  mhc_b                 = "mhc_b",
  mhc_class             = "mhc_class",
  antigen_epitope       = "epitope",
  antigen_gene          = "epitope_gene",
  antigen_species       = "epitope_species",
  reference_id          = "reference_id",
  method_identification = "method_identification",
  method_singlecell     = "method_singlecell",
  method_sequencing     = "method_sequencing",
  method_verification   = "method_verification",
  method_pairing        = "method_pairing",
  meta_study_id         = "study_id",
  meta_cell_subset      = "cell_subset",
  meta_subject_cohort   = "subject_cohort",
  meta_subject_id       = "subject_id",
  meta_replica_id       = "replica_id",
  meta_clone_id         = "clone_id",
  meta_tissue           = "tissue",
  meta_donor_mhc        = "donor_mhc",
  meta_structure_id     = "structure_id"
)

#' Normalise a raw chunk field: strip quoting, then map blanks to NA
#'
#' Two problems are handled together, because they compound.
#'
#' Chunks are read with `quote = ""` so that unescaped quote characters in
#' free-text fields (`comment`) cannot derail the parse. The cost is that
#' genuinely quoted fields arrive with their quote characters intact, so an
#' empty field written as `""` reads as the two-character string `""` rather
#' than as empty. Surrounding quotes are therefore stripped here.
#'
#' Then blanks are mapped to `NA`. Chunks come from hundreds of submitters who
#' spell "missing" as empty, whitespace, or a literal `NA` / `N/A` / `null`.
#' This matters well beyond tidiness: a `meta_subject_id` that survives as a
#' non-missing string inflates donor counts, and donor counts drive every
#' publicness metric in the project.
#'
#' @param x Character vector (coerced).
#' @return Character vector, quote-stripped, with missing spellings as `NA`.
blank_to_na <- function(x) {
  x <- as.character(x)
  ## Strip one layer of matched surrounding quotes, single or double.
  x <- sub('^"(.*)"$', "\\1", x)
  x <- sub("^'(.*)'$", "\\1", x)
  trimmed <- trimws(x)
  x[!is.na(x) & (!nzchar(trimmed) | trimmed %in% c("NA", "na", "N/A", "n/a",
                                                   "NULL", "null", "."))] <- NA_character_
  x
}

#' Locate the VDJdb chunks directory
#'
#' @param dir Directory to search (default `data/raw`).
#' @param negative Return `chunks_negative/` instead of `chunks/`.
#' @return Character scalar path, or `character(0)`.
find_chunks_dir <- function(dir = dir_raw(), negative = FALSE) {
  target <- if (negative) "chunks_negative" else "chunks"
  hits <- list.dirs(dir, recursive = TRUE, full.names = TRUE)
  hits <- hits[basename(hits) == target]
  if (!length(hits)) return(character(0))
  ## Prefer the directory with the most TSVs, in case of stray matches.
  counts <- vapply(hits, function(h) {
    length(list.files(h, pattern = "\\.(tsv|txt)$", ignore.case = TRUE))
  }, integer(1))
  hits[which.max(counts)]
}

#' Read and concatenate VDJdb chunk files
#'
#' Columns vary between chunks (older submissions lack newer fields), so rows
#' are bound with `fill = TRUE` and everything is read as character to stop
#' `fread` inferring incompatible types across files.
#'
#' @param chunks_dir Directory of chunk TSVs.
#' @param verbose Log per-file progress.
#' @return A `data.table` of raw chunk rows, with `chunk_file` added.
read_chunk_files <- function(chunks_dir, verbose = FALSE) {
  files <- list.files(chunks_dir, pattern = "\\.(tsv|txt)$",
                      full.names = TRUE, ignore.case = TRUE)
  files <- files[file.info(files)$size > 0]
  if (!length(files)) stop("No chunk files in ", chunks_dir, call. = FALSE)

  log_msg("reading %d chunk files from %s", length(files), basename(chunks_dir))

  parts <- lapply(files, function(f) {
    dt <- tryCatch(
      data.table::fread(f, sep = "\t", quote = "", colClasses = "character",
                        na.strings = c("", "NA"), showProgress = FALSE),
      error = function(e) {
        log_msg("skipping %s: %s", basename(f), conditionMessage(e), level = "warn")
        NULL
      }
    )
    if (is.null(dt) || !nrow(dt)) return(NULL)
    data.table::setnames(dt, normalise_names(names(dt)))
    if (!"antigen_epitope" %in% names(dt)) return(NULL)

    ## Normalise missingness at the read boundary, for every column, so that
    ## nothing downstream depends on how `fread` happened to treat a blank
    ## field. Chunks are written by many different submitters and spell
    ## "missing" as "", whitespace, or a literal "NA".
    for (cn in names(dt)) {
      data.table::set(dt, j = cn, value = blank_to_na(dt[[cn]]))
    }

    dt[, chunk_file := basename(f)]
    dt
  })

  out <- data.table::rbindlist(Filter(Negate(is.null), parts),
                               fill = TRUE, use.names = TRUE)
  log_msg("concatenated %d rows from %d usable chunks",
          nrow(out), length(Filter(Negate(is.null), parts)))
  out[]
}

#' Load VDJdb from chunk files, in long (one row per chain) form
#'
#' @param chunks_dir Chunks directory; found automatically if `NULL`.
#' @param species Keep only this species. `NULL` keeps all.
#' @param chain "TRB", "TRA", or `NULL` for both.
#' @param negative Read `chunks_negative/` (verified non-binders) instead.
#' @return A `data.table` in the same canonical schema `load_vdjdb()` returns.
load_vdjdb_chunks <- function(chunks_dir = NULL,
                              species = "HomoSapiens",
                              chain = "TRB",
                              negative = FALSE) {

  if (is.null(chunks_dir)) chunks_dir <- find_chunks_dir(negative = negative)
  if (!length(chunks_dir)) {
    stop("No VDJdb chunks directory found under data/raw. ",
         "Run scripts/01_download_data.R, or export from vdjdb.cdr3.net ",
         "(see README).", call. = FALSE)
  }

  raw <- read_chunk_files(chunks_dir)

  ## A stable key per source row, so alpha/beta from one clone stay linked
  ## after reshaping to long form.
  raw[, complex_id := paste(chunk_file, .I, sep = ":")]

  ## --- reshape: one row per chain -------------------------------------------
  ## Columns are assembled into a plain list and converted once, rather than
  ## added with `[[<-` or `:=`, which would copy the table on every column.
  get_col <- function(dt, nm) {
    if (!is.null(nm) && nm %in% names(dt)) blank_to_na(dt[[nm]]) else NA_character_
  }

  build_chain <- function(dt, gene, cdr3_col, v_col, j_col, d_col = NULL) {
    if (!cdr3_col %in% names(dt)) return(NULL)
    cdr3_vals <- blank_to_na(dt[[cdr3_col]])
    keep <- !is.na(cdr3_vals)
    if (!any(keep)) return(NULL)
    sel <- dt[keep]

    cols <- list(
      complex_id = sel$complex_id,
      gene       = gene,
      cdr3       = cdr3_vals[keep],
      v_segm     = get_col(sel, v_col),
      j_segm     = get_col(sel, j_col),
      d_segm     = get_col(sel, d_col),
      chunk_file = sel$chunk_file
    )
    for (src in names(.CHUNK_MAP)) {
      cols[[.CHUNK_MAP[[src]]]] <- get_col(sel, src)
    }
    data.table::as.data.table(cols)
  }

  pieces <- list(
    build_chain(raw, "TRB", "cdr3_beta",  "v_beta",  "j_beta", "d_beta"),
    build_chain(raw, "TRA", "cdr3_alpha", "v_alpha", "j_alpha")
  )
  long <- data.table::rbindlist(Filter(Negate(is.null), pieces),
                                fill = TRUE, use.names = TRUE)

  ## --- donor key -------------------------------------------------------------
  ## `subject_id` is unique only within a study, so the usable donor key is
  ## study (or reference) plus subject.
  long[, has_subject_id := !is.na(subject_id)]
  long[, donor_key := data.table::fifelse(
    has_subject_id,
    paste(data.table::fcoalesce(study_id, reference_id, chunk_file), subject_id, sep = "|"),
    NA_character_
  )]

  ## Chunks carry no assembled score; keep the column for schema compatibility.
  long[, vdjdb_score := NA_integer_]
  long[, is_negative := negative]

  ## --- filters ---------------------------------------------------------------
  n0 <- nrow(long)
  keep_species <- species
  keep_chain   <- chain
  if (!is.null(keep_species)) long <- long[species == keep_species]
  if (!is.null(keep_chain))   long <- long[gene == keep_chain]

  long <- long[!is.na(epitope) & nzchar(epitope)]
  long <- long[grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", cdr3)]

  log_msg("after filtering: %d rows (%.1f%% of %d)",
          nrow(long), 100 * nrow(long) / max(n0, 1L), n0)
  log_msg("  epitopes: %d | donors with IDs: %d | studies: %d",
          data.table::uniqueN(long$epitope),
          data.table::uniqueN(long$donor_key, na.rm = TRUE),
          data.table::uniqueN(long$reference_id))

  long[]
}

#' Summarise how specificities were established
#'
#' Replaces `vdjdb.score` as the basis for confidence filtering. Print this
#' before deciding what to keep - the assay mix differs a lot between studies,
#' and multimer-sorted records are not equivalent to culture-based ones.
#'
#' @param dt Output of `load_vdjdb_chunks()`.
#' @return A `data.table` of method counts, descending.
summarise_methods <- function(dt) {
  if (!"method_identification" %in% names(dt)) {
    log_msg("no method_identification column", level = "warn")
    return(data.table::data.table())
  }
  dt[, .(n_tcr = .N,
         n_epitopes = data.table::uniqueN(epitope)),
     by = .(method_identification)][order(-n_tcr)]
}

#' Flag records whose specificity rests on direct binding evidence
#'
#' Multimer/tetramer sorting and structural data are direct physical evidence
#' of binding. Culture- and activation-based assays are weaker, because a
#' responding cell is not necessarily one that binds the tested peptide.
#' Returns a logical vector rather than filtering, so the caller decides.
#'
#' @param dt Output of `load_vdjdb_chunks()`.
#' @return Logical vector, one element per row.
chunk_confidence <- function(dt) {
  if (!"method_identification" %in% names(dt)) return(rep(NA, nrow(dt)))
  m <- tolower(data.table::fcoalesce(dt$method_identification, ""))
  grepl("multimer|tetramer|dextramer|pentamer|structur|crystal", m)
}

#' Flag TCRs that came from a natural repertoire
#'
#' Some VDJdb records come from in vitro display and selection systems (phage
#' display, yeast display, synthetic libraries) rather than from T cells in a
#' person. Those sequences say nothing about how often recombination produces a
#' given receptor, or how many donors carry it, so they must be excluded from
#' any repertoire-availability or publicness analysis. Because a single display
#' experiment can contribute tens of thousands of sequences against one epitope,
#' leaving them in would dominate exactly the statistics this project measures.
#'
#' They remain legitimate for sequence-motif and structural work, so this
#' returns a flag rather than filtering.
#'
#' @param dt Output of `load_vdjdb_chunks()`.
#' @return Logical vector: TRUE where the TCR plausibly came from a donor.
is_natural_repertoire <- function(dt) {
  if (!"method_identification" %in% names(dt)) return(rep(TRUE, nrow(dt)))
  m <- tolower(data.table::fcoalesce(dt$method_identification, ""))
  !grepl("phage|yeast.?display|synthetic.?librar|in.?vitro.?displa|ribosome.?displa", m)
}
