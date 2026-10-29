## R/publicness_null.R -------------------------------------------------------
## Depth-controlled publicness.
##
## THE PROBLEM THIS SOLVES
##
## `public_fraction` counts CDR3s (or clusters) seen in >= 2 donors. With 81
## donors there are 3,240 donor pairs a sequence could be shared across; with 3
## donors there are 3. Publicness therefore rises with donor count for purely
## combinatorial reasons.
##
## Correlating raw publicness against `n_donors` is close to correlating a
## variable with itself, and in this dataset it produces rho ~ +0.65 whether or
## not any biology is present. Two independent controls are implemented here:
##
##  1. DONOR RAREFACTION - recompute publicness using a fixed number of donors,
##     sampled repeatedly. Every epitope is then measured at the same depth, so
##     values are directly comparable.
##
##  2. PERMUTATION NULL - shuffle donor labels within an epitope. This holds the
##     number of donors, the number of sequences, and each donor's record count
##     exactly fixed, and asks whether similar sequences come from *different*
##     donors more often than chance. The z-score is a depth-free effect size.
##
## The two answer different questions and should agree. Rarefaction makes
## epitopes comparable to each other; the null asks whether any one epitope's
## repertoire is more convergent than its own sampling structure predicts.
## -----------------------------------------------------------------------------

#' Publicness at a fixed number of donors
#'
#' @param cdr3 Character vector of CDR3 sequences.
#' @param donor Character vector of donor keys, same length.
#' @param k Number of donors to sample.
#' @param n_rep Replicates to average over.
#' @param cluster Use cluster-based publicness rather than exact matching.
#' @param max_dist Edit distance for clustering.
#' @param seed Seed, for reproducibility.
#' @return A list with the mean and SD across replicates, or `NA` if the epitope
#'   has fewer than `k` donors - deliberately not extrapolated.
publicness_donor_rarefied <- function(cdr3, donor, k = 5L, n_rep = 50L,
                                      cluster = TRUE, max_dist = 1L,
                                      seed = 1L) {
  ok <- !is.na(cdr3) & !is.na(donor)
  cdr3 <- cdr3[ok]; donor <- donor[ok]
  donors <- unique(donor)

  if (length(donors) < k) {
    return(list(public_rare_mean = NA_real_, public_rare_sd = NA_real_,
                n_rep_used = 0L))
  }

  set.seed(seed)
  vals <- vapply(seq_len(n_rep), function(i) {
    keep <- donor %in% sample(donors, k)
    if (cluster) {
      cluster_publicness(cdr3[keep], donor[keep], max_dist = max_dist)$cluster_public_fraction
    } else {
      publicness(cdr3[keep], donor[keep])$public_fraction
    }
  }, numeric(1))

  vals <- vals[!is.na(vals)]
  if (!length(vals)) {
    return(list(public_rare_mean = NA_real_, public_rare_sd = NA_real_,
                n_rep_used = 0L))
  }
  list(public_rare_mean = mean(vals),
       public_rare_sd   = stats::sd(vals),
       n_rep_used       = length(vals))
}

#' Permutation null for publicness
#'
#' Shuffles donor labels, holding donor count, sequence count and per-donor
#' record counts fixed. The resulting z-score does not depend on sampling depth,
#' so it is comparable across epitopes of very different sizes.
#'
#' Note: with exact matching, a CDR3 observed only once can never be public
#' under any shuffle, so the null is informative mainly for *cluster*
#' publicness, where a cluster holds several distinct sequences. That is why
#' `cluster = TRUE` is the default.
#'
#' @inheritParams publicness_donor_rarefied
#' @param n_perm Number of permutations.
#' @return A list with observed value, null mean and SD, z-score and empirical p.
publicness_null <- function(cdr3, donor, n_perm = 200L, cluster = TRUE,
                            max_dist = 1L, seed = 1L) {
  na_result <- list(observed = NA_real_, null_mean = NA_real_,
                    null_sd = NA_real_, z = NA_real_, p_emp = NA_real_,
                    within_donor_cluster_fraction = NA_real_)

  ok <- !is.na(cdr3) & !is.na(donor)
  cdr3 <- cdr3[ok]; donor <- donor[ok]
  if (length(unique(donor)) < 2L) return(na_result)

  ## --- DEDUPLICATION -------------------------------------------------------
  ## One row per (CDR3, donor). Without this, a donor who contributed a
  ## sequence 500 times would dominate the permutation: shuffling scatters
  ## those 500 records across donors and manufactures sharing that was never
  ## observed. Publicness asks whether a donor carries a sequence at all, not
  ## how many times it was recorded.
  dd <- unique(data.table::data.table(cdr3 = cdr3, donor = donor))
  if (nrow(dd) < 3L) return(na_result)

  ## --- CLUSTER ONCE --------------------------------------------------------
  ## Cluster membership depends only on the sequences, and permutation changes
  ## only the donor labels, so clustering is done once rather than inside every
  ## permutation. (The previous version recomputed it 200 times per epitope,
  ## which is why this step took minutes.)
  cl <- if (cluster) cdr3_clusters(dd$cdr3, max_dist = max_dist) else
    data.table::data.table(cdr3 = unique(dd$cdr3),
                           cluster = seq_along(unique(dd$cdr3)))
  if (is.null(cl)) return(na_result)

  dd <- merge(dd, cl, by = "cdr3", all.x = TRUE)[!is.na(cluster)]
  if (!nrow(dd)) return(na_result)

  ## Public fraction from a vector of donor labels aligned to dd's rows.
  frac_public <- function(donor_vec) {
    tab <- data.table::data.table(cluster = dd$cluster, donor = donor_vec)
    per <- tab[, .(nd = data.table::uniqueN(donor)), by = cluster]
    mean(per$nd >= 2)
  }

  observed <- frac_public(dd$donor)

  ## Diagnostic: how much of the structure is clonal families inside one donor?
  ## A high value explains a negative z-score without implying anything about
  ## cross-donor publicness.
  per_obs <- dd[, .(nd = data.table::uniqueN(donor), n = .N), by = cluster]
  within_frac <- mean(per_obs[n > 1, nd == 1])
  if (is.nan(within_frac)) within_frac <- NA_real_

  set.seed(seed)
  null <- vapply(seq_len(n_perm), function(i) frac_public(sample(dd$donor)),
                 numeric(1))
  null <- null[!is.na(null)]
  if (length(null) < 10L) return(na_result)

  sd_null <- stats::sd(null)
  list(
    observed  = observed,
    null_mean = mean(null),
    null_sd   = sd_null,
    z         = if (sd_null > 0) (observed - mean(null)) / sd_null else NA_real_,
    ## Two-sided: a repertoire can be MORE donor-structured than chance
    ## (clonal families) as well as more public, and both are informative.
    p_emp     = (sum(abs(null - mean(null)) >= abs(observed - mean(null))) + 1) /
                (length(null) + 1),
    within_donor_cluster_fraction = within_frac
  )
}

#' Add depth-controlled publicness to a per-epitope table
#'
#' @param feats Per-epitope features.
#' @param dt TCR-level table.
#' @param k Donors to rarefy to.
#' @param n_rep Rarefaction replicates.
#' @param n_perm Permutations for the null.
#' @param max_dist Clustering distance.
#' @return `feats` with depth-controlled columns joined on `epitope`.
add_depth_controlled <- function(feats, dt, k = 5L, n_rep = 50L,
                                 n_perm = 200L, max_dist = 1L) {
  eps <- feats$epitope
  log_msg("depth-controlled publicness for %d epitopes (k = %d donors, %d perms)",
          length(eps), k, n_perm)

  rows <- lapply(seq_along(eps), function(i) {
    sub <- dt[epitope == eps[i]]
    rare <- publicness_donor_rarefied(sub$cdr3, sub$donor_key, k = k,
                                      n_rep = n_rep, max_dist = max_dist)
    nul  <- publicness_null(sub$cdr3, sub$donor_key, n_perm = n_perm,
                            max_dist = max_dist)
    data.table::data.table(
      epitope            = eps[i],
      public_rare_mean   = rare$public_rare_mean,
      public_rare_sd     = rare$public_rare_sd,
      null_observed      = nul$observed,
      null_mean          = nul$null_mean,
      public_z           = nul$z,
      public_p_emp       = nul$p_emp,
      within_donor_frac  = nul$within_donor_cluster_fraction
    )
  })

  res <- data.table::rbindlist(rows, fill = TRUE)
  out <- merge(feats, res, by = "epitope", all.x = TRUE, sort = FALSE)

  log_msg("  rarefied value available for %d epitopes; z-score for %d",
          sum(!is.na(out$public_rare_mean)), sum(!is.na(out$public_z)))
  out[order(-n_donors, -n_tcr)]
}
