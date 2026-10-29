## scripts/06_figures.R ------------------------------------------------------
## The three figures the README and report carry.
##
## Run:  Rscript scripts/06_figures.R
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

source(file.path("R", "utils.R"))

log_msg("=== 06_figures ===")
ensure_dir(dir_figures())

rp <- function(f) {
  p <- dir_processed(f)
  if (!file.exists(p)) stop(f, " not found - run the earlier scripts first.",
                            call. = FALSE)
  data.table::fread(p, sep = "\t", na.strings = c("", "NA"))
}

depth  <- rp("features_sars2_depth_controlled.tsv")
feat_h <- rp("features_hsv2.tsv")

base_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey35", size = 9.5),
    plot.caption  = element_text(colour = "grey45", size = 8, hjust = 0),
    panel.grid.minor = element_blank()
  )

## --- Figure 1: exact matching misses the sharing that exists -------------------
f1_dt <- data.table::melt(
  depth[, .(epitope, n_donors, protein,
            `exact match` = public_fraction,
            `cluster, dist <= 1` = cluster_public_fraction,
            `cluster, dist <= 2` = cluster_public_fraction_d2)],
  id.vars = c("epitope", "n_donors", "protein"),
  variable.name = "measure", value.name = "publicness"
)[!is.na(publicness)]

p1 <- ggplot(f1_dt, aes(measure, publicness)) +
  geom_boxplot(outlier.shape = NA, fill = "grey93", colour = "grey55", width = 0.5) +
  geom_jitter(aes(colour = protein), width = 0.12, height = 0, size = 2.2, alpha = 0.85) +
  labs(
    title    = "Exact CDR3 matching understates repertoire sharing",
    subtitle = "People responding to one epitope converge on motifs, not identical sequences",
    x = NULL, y = "Public fraction", colour = "Protein",
    caption = "Each point is one SARS-CoV-2 structural epitope (n = 17)."
  ) +
  base_theme
ggsave(dir_figures("fig1_exact_vs_cluster.png"), p1, width = 7.5, height = 4.6, dpi = 150)

## --- Figure 2: the depth confound, and what survives it ------------------------
f2_dt <- data.table::melt(
  depth[, .(epitope, n_donors, protein,
            `uncontrolled` = cluster_public_fraction,
            `rarefied to 3 donors` = public_rare_mean)],
  id.vars = c("epitope", "n_donors", "protein"),
  variable.name = "measure", value.name = "publicness"
)[!is.na(publicness)]

labs_dt <- f2_dt[, .(lab = sprintf("rho = %+.2f",
                                   stats::cor(publicness, n_donors,
                                              method = "spearman"))),
                 by = measure]

p2 <- ggplot(f2_dt, aes(n_donors, publicness)) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey70",
              linewidth = 0.6, formula = y ~ x) +
  geom_point(aes(colour = protein), size = 2.4, alpha = 0.9) +
  geom_text(data = labs_dt, aes(label = lab), x = Inf, y = Inf,
            hjust = 1.1, vjust = 1.6, size = 3.4, colour = "grey25") +
  facet_wrap(~measure) +
  scale_x_log10() +
  labs(
    title    = "The association survives depth control",
    subtitle = "Raw publicness rises with donor count combinatorially; rarefying to a fixed donor count removes that",
    x = "Donors (log scale)", y = "Public fraction", colour = "Protein",
    caption = "Spearman correlations, n = 17. No p-values: the sample is too small to support them."
  ) +
  base_theme
ggsave(dir_figures("fig2_depth_control.png"), p2, width = 9, height = 4.4, dpi = 150)

## --- Figure 3: the data-scarcity finding ---------------------------------------
f3 <- rbind(
  depth[, .(label = epitope, n_donors, n_tcr, system = "SARS-CoV-2 structural")],
  feat_h[, .(label = epitope,
             n_donors = data.table::fifelse(is.na(n_donors), 0L, n_donors),
             n_tcr,
             system = data.table::fifelse(protein == "gD",
                                          "HSV-2 glycoprotein D",
                                          "HSV-2 (other proteins)"))]
)
f3[, system := factor(system, levels = c("SARS-CoV-2 structural",
                                         "HSV-2 (other proteins)",
                                         "HSV-2 glycoprotein D"))]

p3 <- ggplot(f3, aes(n_tcr, pmax(n_donors, 0.5), colour = system)) +
  geom_point(size = 2.8, alpha = 0.85) +
  ggrepel::geom_text_repel(
    data = f3[system == "HSV-2 glycoprotein D"],
    aes(label = label), size = 3.2, nudge_y = 0.4, show.legend = FALSE
  ) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c("SARS-CoV-2 structural" = "grey70",
                                 "HSV-2 (other proteins)" = "#4C72B0",
                                 "HSV-2 glycoprotein D" = "#C44E52")) +
  labs(
    title    = "Repertoire data exists in proportion to attention, not clinical need",
    subtitle = "Glycoprotein D, a leading genital herpes vaccine target, has one epitope and four TCRs in VDJdb",
    x = "TCR records (log scale)", y = "Donors (log scale, 0 shown at 0.5)",
    colour = NULL,
    caption = "No repertoire-availability metric is computable at the gD sample size."
  ) +
  base_theme + theme(legend.position = "bottom")

if (requireNamespace("ggrepel", quietly = TRUE)) {
  ggsave(dir_figures("fig3_data_scarcity.png"), p3, width = 8, height = 5.2, dpi = 150)
} else {
  log_msg("ggrepel not installed; drawing fig3 without labels", level = "warn")
  p3b <- p3; p3b$layers <- p3b$layers[1]
  ggsave(dir_figures("fig3_data_scarcity.png"), p3b, width = 8, height = 5.2, dpi = 150)
}

log_msg("wrote 3 figures to figures/", level = "ok")
log_msg("done.", level = "ok")
