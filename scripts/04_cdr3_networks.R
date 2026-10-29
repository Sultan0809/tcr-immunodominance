## scripts/04_cdr3_networks.R -------------------------------------------------
## Cluster-based publicness, and how it compares to exact-match publicness.
##
## Run:  Rscript scripts/04_cdr3_networks.R
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path("R", "utils.R"))
source(file.path("R", "diversity.R"))
source(file.path("R", "features.R"))
source(file.path("R", "cdr3_network.R"))

log_msg("=== 04_cdr3_networks ===")

read_processed <- function(f) {
  p <- dir_processed(f)
  if (!file.exists(p)) stop(f, " not found - run the earlier scripts first.",
                            call. = FALSE)
  data.table::fread(p, sep = "\t", na.strings = c("", "NA"))
}

sars2      <- read_processed("sars2_trb.tsv")
analysable <- read_processed("features_sars2_analysable.tsv")

if (!nrow(analysable)) stop("No analysable epitopes.", call. = FALSE)

## --- 1. Cluster publicness at two thresholds -----------------------------------
## max_dist = 1 is conservative; 2 is common in the literature and much more
## permissive. Reporting both shows how sensitive the conclusion is to a
## parameter that is usually chosen by convention rather than by evidence.
res1 <- add_cluster_features(analysable, sars2, max_dist = 1L)
res2 <- add_cluster_features(analysable, sars2, max_dist = 2L)

data.table::setnames(res2,
  c("cluster_public_fraction", "mean_donors_per_cluster", "n_clusters",
    "largest_cluster_fraction"),
  c("cluster_public_fraction_d2", "mean_donors_per_cluster_d2",
    "n_clusters_d2", "largest_cluster_fraction_d2"))

combined <- merge(
  res1,
  res2[, .(epitope, cluster_public_fraction_d2, mean_donors_per_cluster_d2,
           n_clusters_d2, largest_cluster_fraction_d2)],
  by = "epitope", sort = FALSE
)
write_processed(combined, "features_sars2_with_clusters.tsv")

## --- 2. Does clustering recover publicness that exact matching misses? ----------
log_msg("--- exact vs cluster publicness ---")
cmp <- combined[, .(
  epitope, protein, n_donors, n_cdr3,
  exact      = round(public_fraction, 3),
  clust_d1   = round(cluster_public_fraction, 3),
  clust_d2   = round(cluster_public_fraction_d2, 3),
  n_clusters, n_clusters_d2,
  chain_d2   = round(largest_cluster_fraction_d2, 3)
)][order(-n_donors)]
print(cmp)

message("\nMedian public fraction:")
message("  exact match:        ", round(median(combined$public_fraction, na.rm = TRUE), 4))
message("  cluster, dist <= 1: ", round(median(combined$cluster_public_fraction, na.rm = TRUE), 4))
message("  cluster, dist <= 2: ", round(median(combined$cluster_public_fraction_d2, na.rm = TRUE), 4))

## Chaining check: a component that swallows most sequences inflates publicness.
bad <- combined[largest_cluster_fraction_d2 > 0.5, .N]
if (bad > 0) {
  log_msg("%d epitope(s) have one cluster holding >50%% of sequences at dist<=2;",
          bad, level = "warn")
  log_msg("  single-linkage chaining - prefer dist<=1 for those.", level = "warn")
}

## --- 3. Is repertoire availability associated with donor breadth? ---------------
## The project's central question, on the data that exist. With ~17 epitopes
## this is descriptive: Spearman on a sample this size is unstable, so the
## correlation is reported with its n and without a p-value.
log_msg("--- association with donor breadth (descriptive only) ---")
for (v in c("public_fraction", "cluster_public_fraction",
            "cluster_public_fraction_d2", "v_entropy", "clonality")) {
  x <- combined[[v]]
  ok <- !is.na(x) & !is.na(combined$n_donors)
  if (sum(ok) >= 5) {
    rho <- stats::cor(x[ok], combined$n_donors[ok], method = "spearman")
    message(sprintf("  %-28s rho = %+.3f  (n = %d)", v, rho, sum(ok)))
  } else {
    message(sprintf("  %-28s too few epitopes (n = %d)", v, sum(ok)))
  }
}

message("\nNOTE: n is small and these epitopes differ in HLA restriction, which")
message("is itself a determinant of dominance. Treat as exploratory, not a test.")

log_msg("done.", level = "ok")
