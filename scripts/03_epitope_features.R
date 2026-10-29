## scripts/03_epitope_features.R ----------------------------------------------
## Build the per-epitope feature table and report what is analysable.
##
## Run:  Rscript scripts/03_epitope_features.R
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path("R", "utils.R"))
source(file.path("R", "diversity.R"))
source(file.path("R", "features.R"))

log_msg("=== 03_epitope_features ===")

read_processed <- function(f) {
  p <- dir_processed(f)
  if (!file.exists(p)) {
    stop(f, " not found - run scripts/02_parse_vdjdb.R first.", call. = FALSE)
  }
  data.table::fread(p, sep = "\t", na.strings = c("", "NA"))
}

sars2 <- read_processed("sars2_trb.tsv")
hsv2  <- read_processed("hsv2_trb.tsv")

## --- 1. Features ---------------------------------------------------------------
## rarefy_depth is a trade-off: deeper is a better comparison but keeps fewer
## epitopes. 30 is a starting point - the sensitivity check below shows the cost.
RAREFY_DEPTH <- 30L

log_msg("--- SARS-CoV-2 ---")
feat_sars2 <- epitope_features(sars2, rarefy_depth = RAREFY_DEPTH, min_tcr = 10L)

log_msg("--- HSV-2 ---")
feat_hsv2 <- epitope_features(hsv2, rarefy_depth = RAREFY_DEPTH, min_tcr = 1L)

write_processed(feat_sars2, "features_sars2.tsv")
write_processed(feat_hsv2,  "features_hsv2.tsv")

## --- 2. Analysable subset ------------------------------------------------------
log_msg("--- analysable subset (structural proteins) ---")
struct <- feat_sars2[structural %in% TRUE]
analysable <- analysable_epitopes(struct, min_donors = 3L, min_tcr = 20L)
write_processed(analysable, "features_sars2_analysable.tsv")

## --- 3. How much does the depth threshold cost? ---------------------------------
## Rarefaction depth trades comparability against sample size. Make that
## trade-off explicit rather than hiding it in a default.
log_msg("--- sensitivity of the analysable set to thresholds ---")
## Features are recomputed once per depth, not once per cell.
depths <- c(10L, 20L, 30L, 50L, 100L)
by_depth <- lapply(depths, function(d) {
  suppressMessages(
    epitope_features(sars2[structural %in% TRUE], rarefy_depth = d, min_tcr = 10L)
  )
})
names(by_depth) <- as.character(depths)

grid <- data.table::CJ(depth = depths, min_donors = c(2L, 3L, 5L))
grid[, n_epitopes := {
  f <- by_depth[[as.character(depth)]]
  nrow(f[!is.na(n_donors) & n_donors >= min_donors & !is.na(richness_rare)])
}, by = .(depth, min_donors)]
print(dcast(grid, depth ~ min_donors, value.var = "n_epitopes"))
message("(rows = rarefaction depth, columns = minimum donors, cells = epitopes retained)")

## --- 4. Report ------------------------------------------------------------------
log_msg("--- top epitopes by donor breadth ---")
cols <- c("epitope", "protein", "mhc_a", "n_tcr", "n_cdr3", "n_donors",
          "n_studies", "public_fraction", "mean_donors_per_cdr3",
          "richness_rare", "v_entropy", "clonality")
cols <- intersect(cols, names(analysable))
if (nrow(analysable)) print(head(analysable[, ..cols], 20))

## The core contrast of the project, with the record-count confound removed.
if (nrow(analysable) >= 4) {
  log_msg("--- S vs N (analysable epitopes only) ---")
  ## as.numeric() throughout: median() of an integer vector returns integer for
  ## odd-length groups and double for even, and data.table requires one type
  ## across groups.
  print(analysable[protein %in% c("S", "N"), .(
    n_epitopes       = .N,
    med_donors       = as.numeric(median(n_donors)),
    med_public_frac  = round(as.numeric(median(public_fraction, na.rm = TRUE)), 3),
    med_rare_rich    = round(as.numeric(median(richness_rare, na.rm = TRUE)), 1),
    med_v_entropy    = round(as.numeric(median(v_entropy, na.rm = TRUE)), 3),
    med_clonality    = round(as.numeric(median(clonality, na.rm = TRUE)), 3)
  ), by = protein])
}

log_msg("--- HSV-2 ---")
hcols <- intersect(cols, names(feat_hsv2))
if (nrow(feat_hsv2)) print(feat_hsv2[, ..hcols])

gd <- feat_hsv2[protein == "gD"]
message("\ngD epitopes in VDJdb: ", nrow(gd))
if (nrow(gd)) {
  message("  total TCRs: ", sum(gd$n_tcr),
          " | donors: ", sum(gd$n_donors, na.rm = TRUE),
          " | studies: ", sum(gd$n_studies))
  message("  No repertoire-availability metric is computable at this scale.")
  message("  This is the step-4 finding, not a pipeline failure.")
}

log_msg("done.", level = "ok")
