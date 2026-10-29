## scripts/01_download_data.R -------------------------------------------------
## Acquire raw inputs: VDJdb + UniProt reference proteins.
## Idempotent - safe to re-run; nothing is re-downloaded if already cached.
##
## Run:  Rscript scripts/01_download_data.R
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path("R", "utils.R"))
source(file.path("R", "load_vdjdb.R"))
source(file.path("R", "annotate_antigen.R"))

log_msg("=== 01_download_data ===")

ensure_dir(dir_raw())
ensure_dir(dir_processed())

## --- 1. VDJdb ----------------------------------------------------------------
got_vdjdb <- download_vdjdb()

if (!got_vdjdb) {
  log_msg(paste(
    "Continuing without VDJdb. Script 02 will stop until a file is present.",
    "See the README section 'If the VDJdb auto-download fails'.",
    sep = "\n         "
  ), level = "warn")
}

## --- 2. Reference proteins ---------------------------------------------------
## Fetched individually and cached as FASTA under data/raw/uniprot/.
proteome <- build_reference_proteome()

write_processed(
  proteome[, .(organism, protein, structural, accession, length, header)],
  "reference_proteome_index.tsv"
)

## The sequences themselves are kept separately - they are the only large
## column and nothing downstream needs them joined to the index.
saveRDS(proteome, dir_processed("reference_proteome.rds"))
log_msg("wrote reference_proteome.rds", level = "ok")

## --- summary -----------------------------------------------------------------
log_msg("--- inputs ---")
print(proteome[, .(organism, protein, structural, accession, length)])

log_msg("done. Next: Rscript scripts/02_parse_vdjdb.R", level = "ok")
