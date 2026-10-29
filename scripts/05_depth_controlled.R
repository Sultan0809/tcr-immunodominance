## scripts/05_depth_controlled.R ----------------------------------------------
## Does repertoire availability track donor breadth once sampling depth is
## controlled for? This is the project's central question, tested properly.
##
## Run:  Rscript scripts/05_depth_controlled.R
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

source(file.path("R", "utils.R"))
source(file.path("R", "diversity.R"))
source(file.path("R", "features.R"))
source(file.path("R", "cdr3_network.R"))
source(file.path("R", "publicness_null.R"))

log_msg("=== 05_depth_controlled ===")

read_processed <- function(f) {
  p <- dir_processed(f)
  if (!file.exists(p)) stop(f, " not found - run the earlier scripts first.",
                            call. = FALSE)
  data.table::fread(p, sep = "\t", na.strings = c("", "NA"))
}

sars2 <- read_processed("sars2_trb.tsv")
feats <- read_processed("features_sars2_with_clusters.tsv")

## --- 1. Depth-controlled publicness ---------------------------------------------
## k = 3 keeps every analysable epitope (the smallest has 3 donors); k = 5 is a
## firmer comparison but drops a few. Both are reported.
res <- add_depth_controlled(feats, sars2, k = 3L, n_rep = 50L, n_perm = 200L)

res5 <- add_depth_controlled(feats[n_donors >= 5], sars2, k = 5L,
                             n_rep = 50L, n_perm = 200L)
data.table::setnames(res5, "public_rare_mean", "public_rare_k5")
res <- merge(res, res5[, .(epitope, public_rare_k5)],
             by = "epitope", all.x = TRUE, sort = FALSE)

write_processed(res, "features_sars2_depth_controlled.tsv")

## --- 2. The comparison that matters ----------------------------------------------
log_msg("--- raw vs depth-controlled publicness ---")
tbl <- res[, .(
  epitope, protein, n_donors,
  raw_cluster = round(cluster_public_fraction, 3),
  rare_k3     = round(public_rare_mean, 3),
  rare_k5     = round(public_rare_k5, 3),
  z           = round(public_z, 2),
  p           = round(public_p_emp, 3),
  within_don  = round(within_donor_frac, 2)
)][order(-n_donors)]
print(tbl)

message("\n`within_don` = fraction of multi-sequence clusters confined to ONE donor.")
message("High values mean within-donor clonal families dominate the similarity")
message("structure, which drives z negative independently of cross-donor sharing.")

## --- 3. Correlations, before and after control ------------------------------------
log_msg("--- Spearman vs n_donors ---")
cor_vs_donors <- function(v, label) {
  x <- res[[v]]
  ok <- !is.na(x) & !is.na(res$n_donors)
  if (sum(ok) < 5) {
    message(sprintf("  %-26s too few epitopes (n = %d)", label, sum(ok)))
    return(invisible(NULL))
  }
  rho <- stats::cor(x[ok], res$n_donors[ok], method = "spearman")
  message(sprintf("  %-26s rho = %+.3f  (n = %d)", label, rho, sum(ok)))
}

message("\nUncontrolled (expected to be inflated by sampling depth):")
cor_vs_donors("public_fraction", "exact publicness")
cor_vs_donors("cluster_public_fraction", "cluster publicness")

message("\nDepth-controlled (the honest estimates):")
cor_vs_donors("public_rare_mean", "rarefied to 3 donors")
cor_vs_donors("public_rare_k5", "rarefied to 5 donors")
cor_vs_donors("public_z", "permutation z-score")

message("\nINTERPRETATION")
message("  If the uncontrolled correlations are strong and the controlled ones")
message("  collapse toward zero, the apparent association was sampling depth,")
message("  not biology - a negative result, and a real one worth reporting.")
message("  If the controlled estimates hold up, repertoire convergence tracks")
message("  donor breadth independently of how much each epitope was studied.")

## --- 4. Is any epitope's repertoire more convergent than chance? -------------------
log_msg("--- repertoires more DONOR-STRUCTURED than chance (z < -2) ---")
ds <- res[!is.na(public_z) & public_z < -2][order(public_z)]
if (nrow(ds)) {
  message("  Similar sequences cluster within donors more than expected:")
  print(ds[, .(epitope, protein, n_donors,
               observed = round(null_observed, 3),
               expected = round(null_mean, 3),
               z = round(public_z, 2),
               within_don = round(within_donor_frac, 2))])
  message("  This is clonal expansion inside donors, not evidence about publicness.")
} else {
  message("  None.")
}

log_msg("--- epitopes with convergent repertoires (z > 2) ---")
conv <- res[!is.na(public_z) & public_z > 2][order(-public_z)]
if (nrow(conv)) {
  print(conv[, .(epitope, protein, mhc_a, n_donors,
                 observed = round(null_observed, 3),
                 expected = round(null_mean, 3),
                 z = round(public_z, 2), p = round(public_p_emp, 3))])
} else {
  message("  None. At this sample depth, cross-donor sharing is consistent with")
  message("  the sampling structure alone.")
}

## --- 5. Figure ---------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  ensure_dir(dir_figures())

  plot_dt <- data.table::melt(
    res[, .(epitope, n_donors, protein,
            `raw cluster` = cluster_public_fraction,
            `rarefied (k=3)` = public_rare_mean)],
    id.vars = c("epitope", "n_donors", "protein"),
    variable.name = "measure", value.name = "publicness"
  )[!is.na(publicness)]

  p <- ggplot(plot_dt, aes(n_donors, publicness, colour = protein)) +
    geom_point(size = 2.4, alpha = 0.85) +
    facet_wrap(~measure, scales = "free_y") +
    scale_x_log10() +
    labs(
      title    = "Publicness against donor breadth, before and after depth control",
      subtitle = "Raw publicness rises with donor count for combinatorial reasons; rarefying to a fixed donor count removes that",
      x        = "Donors (log scale)", y = "Public fraction", colour = "Protein"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.subtitle = element_text(colour = "grey35"))

  ggsave(dir_figures("publicness_depth_control.png"), p,
         width = 9, height = 4.2, dpi = 150)
  log_msg("wrote figures/publicness_depth_control.png", level = "ok")
}

log_msg("done.", level = "ok")
