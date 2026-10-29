## scripts/02_parse_vdjdb.R ---------------------------------------------------
## Parse VDJdb, annotate epitopes against the reference proteome, write the
## tidy tables that steps 2-4 consume, and emit a QC summary.
##
## Prefers VDJdb's per-study chunk files (richer: real meta_* columns, paired
## chains, verified non-binders) and falls back to an assembled table or web
## export if chunks are not present.
##
## Run:  Rscript scripts/02_parse_vdjdb.R
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path("R", "utils.R"))
source(file.path("R", "load_vdjdb.R"))
source(file.path("R", "load_vdjdb_chunks.R"))
source(file.path("R", "annotate_antigen.R"))

log_msg("=== 02_parse_vdjdb ===")

## --- 1. Load ------------------------------------------------------------------
chunks_dir <- find_chunks_dir()

if (length(chunks_dir)) {
  log_msg("source: VDJdb chunk files")
  vdjdb <- load_vdjdb_chunks(chunks_dir, species = "HomoSapiens", chain = "TRB")
  source_kind <- "chunks"
} else {
  log_msg("source: assembled VDJdb table / web export")
  vdjdb <- load_vdjdb(species = "HomoSapiens", chain = "TRB", min_score = 1L)
  source_kind <- "assembled"
}

proteome_path <- dir_processed("reference_proteome.rds")
if (!file.exists(proteome_path)) {
  stop("reference_proteome.rds not found - run scripts/01_download_data.R first.",
       call. = FALSE)
}
proteome <- readRDS(proteome_path)

## --- 2. Confidence -------------------------------------------------------------
## Chunks carry no assembled vdjdb.score, so confidence comes from how the
## specificity was established. Nothing is dropped here - the flag is recorded
## and step 2 decides, after the distribution below has been inspected.
if (source_kind == "chunks") {
  vdjdb[, direct_binding := chunk_confidence(vdjdb)]
  vdjdb[, natural_repertoire := is_natural_repertoire(vdjdb)]

  log_msg("--- how specificities were established ---")
  print(head(summarise_methods(vdjdb), 15))
  log_msg("records with direct binding evidence (multimer/structure): %d (%.1f%%)",
          sum(vdjdb$direct_binding, na.rm = TRUE),
          100 * mean(vdjdb$direct_binding, na.rm = TRUE))

  ## In vitro display sequences are not from any donor's repertoire. A single
  ## such experiment can contribute tens of thousands of sequences against one
  ## epitope, which would dominate the very statistics this project measures.
  n_display <- sum(!vdjdb$natural_repertoire, na.rm = TRUE)
  if (n_display > 0) {
    log_msg("excluding %d in vitro display records (%.1f%%) from repertoire analysis",
            n_display, 100 * mean(!vdjdb$natural_repertoire, na.rm = TRUE),
            level = "warn")
    print(vdjdb[natural_repertoire %in% FALSE,
                .(n_tcr = .N, n_epitopes = uniqueN(epitope)),
                by = method_identification][order(-n_tcr)])
    vdjdb <- vdjdb[natural_repertoire %in% TRUE]
  }
} else {
  vdjdb[, direct_binding := NA]
  vdjdb[, natural_repertoire := TRUE]
}

## --- 3. Annotate ---------------------------------------------------------------
annotated <- annotate_vdjdb(vdjdb, proteome)

## --- 4. Split out the two study systems ---------------------------------------
sars2 <- annotated[
  organism == "SARS-CoV-2" |
    grepl("SARS.?CoV.?2|SARS2", epitope_species, ignore.case = TRUE)
]
hsv2 <- annotated[
  organism == "HSV-2" |
    grepl("herpes|HSV", epitope_species, ignore.case = TRUE)
]
sars2_structural <- sars2[structural %in% TRUE]

## --- 5. Write ------------------------------------------------------------------
keep_cols <- c(
  "epitope", "organism", "protein", "protein_label", "structural",
  "start", "end", "rel_position",
  "cdr3", "v_segm", "j_segm", "d_segm", "gene",
  "mhc_a", "mhc_b", "mhc_class",
  "method_identification", "direct_binding", "natural_repertoire",
  "reference_id", "study_id", "subject_id", "donor_key", "has_subject_id",
  "cell_subset", "tissue", "complex_id"
)
keep_cols <- intersect(keep_cols, names(annotated))

write_processed(annotated[, ..keep_cols],        "vdjdb_human_trb.tsv")
write_processed(sars2[, ..keep_cols],            "sars2_trb.tsv")
write_processed(sars2_structural[, ..keep_cols], "sars2_structural_trb.tsv")
write_processed(hsv2[, ..keep_cols],             "hsv2_trb.tsv")

## --- 6. QC summary -------------------------------------------------------------
qc_by_epitope <- function(dt) {
  if (!nrow(dt)) return(data.table::data.table())
  dt[, {
    mhcs <- sort(unique(mhc_a[!is.na(mhc_a)]))
    .(
      n_tcr         = .N,
      n_unique_cdr3 = data.table::uniqueN(cdr3),
      n_donors      = data.table::uniqueN(donor_key[!is.na(donor_key)]),
      n_studies     = data.table::uniqueN(reference_id),
      pct_direct    = round(100 * mean(direct_binding, na.rm = TRUE), 1),
      mhc           = if (length(mhcs)) paste(utils::head(mhcs, 2L), collapse = "; ") else NA_character_,
      protein       = protein_label[1]
    )
  }, by = epitope][order(-n_tcr)]
}

qc_sars2 <- qc_by_epitope(sars2_structural)
qc_hsv2  <- qc_by_epitope(hsv2)

if (nrow(qc_sars2)) write_processed(qc_sars2, "qc_epitopes_sars2_structural.tsv")
if (nrow(qc_hsv2))  write_processed(qc_hsv2,  "qc_epitopes_hsv2.tsv")

## --- 7. Report -----------------------------------------------------------------
log_msg("--- QC summary ---")
message("\nAll human TRB records:        ", nrow(annotated))
message("SARS-CoV-2 (any protein):     ", nrow(sars2),
        "  across ", uniqueN(sars2$epitope), " epitopes")
message("SARS-CoV-2 structural (SNME): ", nrow(sars2_structural),
        "  across ", uniqueN(sars2_structural$epitope), " epitopes")
message("HSV-2 (any protein):          ", nrow(hsv2),
        "  across ", uniqueN(hsv2$epitope), " epitopes")

if (nrow(sars2_structural)) {
  message("\nRecords per structural protein:")
  print(sars2_structural[, .(n_tcr = .N, n_epitopes = uniqueN(epitope)),
                         by = protein][order(-n_tcr)])
  message("\nTop 15 SARS-CoV-2 structural epitopes by TCR count:")
  print(head(qc_sars2, 15))
}

message("\nHSV-2 availability (the step-4 case study):")
if (nrow(qc_hsv2)) {
  print(qc_hsv2)
} else {
  message("  No HSV-2 TRB records passed filtering. This is the expected ",
          "result and is itself the finding for step 4.")
}

cov <- annotated[, mean(has_subject_id)]
log_msg("donor IDs available for %.1f%% of records", 100 * cov)
log_msg("  publicness metrics in step 2 are restricted to that subset.")

log_msg("done.", level = "ok")
