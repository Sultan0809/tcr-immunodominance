## scripts/07_predicted_epitopes.R --------------------------------------------
## Do the epitopes from an external prediction pipeline have observed human TCR
## repertoires, and how broad are they?
##
## Input:  data/input/epitopes_predicted.csv  (columns: epitope, organism,
##         mhc_class, source_table, position, iedb_match_type)
## Run:    Rscript scripts/07_predicted_epitopes.R
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path("R", "utils.R"))
source(file.path("R", "diversity.R"))
source(file.path("R", "features.R"))
source(file.path("R", "match_epitopes.R"))

log_msg("=== 07_predicted_epitopes ===")

rp <- function(f) {
  p <- dir_processed(f)
  if (!file.exists(p)) stop(f, " not found - run the earlier scripts first.",
                            call. = FALSE)
  data.table::fread(p, sep = "\t", na.strings = c("", "NA"))
}

pred_path <- proj_path("data", "input", "epitopes_predicted.csv")
if (!file.exists(pred_path)) {
  stop("data/input/epitopes_predicted.csv not found.", call. = FALSE)
}
predicted <- data.table::fread(pred_path, na.strings = c("", "NA"))
predicted[, epitope := toupper(trimws(epitope))]
log_msg("read %d predicted epitopes (%d unique)",
        nrow(predicted), uniqueN(predicted$epitope))

## --- 1. Observed repertoires ---------------------------------------------------
## IMPORTANT: match against EVERY epitope that has any sequenced TCR, not just
## the analysable subset. The analysable set is filtered to >= 3 donors and
## >= 20 records because the publicness statistics need that depth - but an
## epitope with 2 TCRs still answers "has a human TCR ever been sequenced
## against this?", which is the question here. Matching against the filtered
## set would report real data as absent.
tcr_level <- rbind(rp("sars2_trb.tsv"), rp("hsv2_trb.tsv"), fill = TRUE)

observed <- epitope_features(tcr_level, rarefy_depth = 30L, min_tcr = 1L)
log_msg("observed set: %d epitopes with at least one sequenced TCR",
        nrow(observed))

## Join the depth-controlled publicness measures where they exist. Most
## epitopes are too shallow for them; those columns stay NA, which is correct.
depth_path <- dir_processed("features_sars2_depth_controlled.tsv")
if (file.exists(depth_path)) {
  dc <- rp("features_sars2_depth_controlled.tsv")
  dc_cols <- intersect(c("epitope", "cluster_public_fraction",
                         "public_rare_mean", "public_z"), names(dc))
  observed <- merge(observed, dc[, ..dc_cols], by = "epitope",
                    all.x = TRUE, sort = FALSE)
  log_msg("  depth-controlled measures attached for %d of them",
          sum(!is.na(observed$public_rare_mean)))
}

## --- 2. Match -------------------------------------------------------------------
res <- match_predicted_epitopes(predicted, observed, min_core = 7L)
write_processed(res, "predicted_epitope_matches.tsv")

## --- 3. Report -------------------------------------------------------------------
show_cols <- intersect(
  c("epitope", "mhc_class", "matched_epitope", "match_type", "n_subs",
    "n_tcr", "n_donors", "cluster_public_fraction", "public_rare_mean"),
  names(res)
)

for (org in unique(res$organism)) {
  sub <- res[organism == org]
  hit <- sub[!is.na(match_type)]

  log_msg("--- %s ---", org)
  message(sprintf("  %d of %d predicted epitopes have an observed TCR repertoire",
                  uniqueN(hit$epitope), uniqueN(sub$epitope)))

  if (nrow(hit)) {
    for (cls in c("same epitope", "point variant", "overlapping")) {
      h <- hit[match_class %in% cls]
      if (!nrow(h)) next
      message(sprintf("\n  [%s] %d epitope(s):", cls, uniqueN(h$epitope)))
      if (cls == "overlapping") {
        ## Weakest evidence: often a different epitope from the same region,
        ## with its own HLA restriction. Listed, never counted as recognition.
        message("    ", paste(sort(unique(h$epitope)), collapse = ", "))
      } else {
        print(h[, ..show_cols])
      }
    }

    ## Separate "a TCR exists" from "enough data to characterise the
    ## repertoire" - they support very different claims.
    deep <- hit[match_class %in% "same epitope" &
                  !is.na(n_donors) & n_donors >= 3 & n_tcr >= 20]
    message(sprintf("\n  Of these, %d have enough depth (>= 3 donors, >= 20 records)",
                    uniqueN(deep$epitope)))
    message("  for the repertoire statistics; the rest confirm recognition only.")
  }

  miss <- setdiff(sub$epitope, hit$epitope)
  if (length(miss)) {
    message("\n  No observed TCRs (", length(miss), "):")
    message("    ", paste(miss, collapse = ", "))
  }
  message("")
}

## --- 4. Where prediction and repertoire breadth agree ------------------------------
log_msg("--- predicted epitopes ranked by observed donor breadth ---")
ranked <- res[match_class %in% "same epitope" &
                !is.na(n_donors) & n_donors > 0][order(-n_donors)]
if (nrow(ranked)) {
  print(unique(ranked[, .(epitope, organism, mhc_class, matched_epitope,
                          match_type, n_donors, n_tcr,
                          publicness = round(cluster_public_fraction, 3))]))
} else {
  message("  None of the predicted epitopes has donor-level repertoire data.")
}

## --- 5. Summary ---------------------------------------------------------------------
log_msg("--- summary ---")
## NB: `%in%` not `==`. For unmatched epitopes match_class is NA, and
## `NA == "same epitope"` yields NA, which selects NA elements that uniqueN()
## then counts as a distinct value - inflating every cell by one.
summ <- res[, .(
  predicted        = uniqueN(epitope),
  same_epitope     = uniqueN(epitope[match_class %in% "same epitope"]),
  point_variant    = uniqueN(epitope[match_class %in% "point variant"]),
  any_overlap      = uniqueN(epitope[!is.na(match_type)]),
  with_donor_ids   = uniqueN(epitope[match_class %in% "same epitope" &
                                       !is.na(n_donors) & n_donors > 0]),
  analysable       = uniqueN(epitope[match_class %in% "same epitope" &
                                       !is.na(n_donors) & n_donors >= 3 & n_tcr >= 20])
), by = .(organism, mhc_class)][order(organism, mhc_class)]

## Internal consistency: a class can never exceed the total matched.
stopifnot(all(summ$same_epitope <= summ$any_overlap),
          all(summ$analysable <= summ$with_donor_ids))
print(summ)

message("\nA miss means no sequenced human TCR exists for that epitope.")
message("It is NOT evidence the epitope is non-immunogenic - IEDB assay")
message("validation and TCR sequencing measure different things.")

log_msg("done.", level = "ok")
