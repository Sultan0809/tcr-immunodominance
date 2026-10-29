## R/features.R --------------------------------------------------------------
## Per-epitope repertoire features: the core table of the project.
##
## Design decisions forced by the data (see reports/):
##
##  * `n_donors` is the dominance proxy, NOT `n_tcr`. LLLDRLNQL has 1753 TCRs
##    from 3 donors; TFEYVSQPFLMDLE has 2220 from 70. Record counts measure
##    study effort; donor counts measure how broadly a response was observed.
##
##  * VDJdb chunks carry amino-acid CDR3s only, with no nucleotide column, so
##    convergent recombination cannot be computed. Cross-donor publicness is
##    therefore the sole repertoire-availability proxy available in R, which is
##    why OLGA generation probability is worth adding when possible.
##
##  * Every diversity measure is confounded by sample size, so each is reported
##    both raw and rarefied to a common depth.
## -----------------------------------------------------------------------------

#' Normalised Shannon entropy of a categorical vector
#'
#' Used for V and J gene usage: 0 means one gene dominates, 1 means uniform.
#'
#' @param x Character vector (NAs dropped).
#' @return Numeric in [0, 1], or `NA` if empty.
usage_entropy <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_real_)
  shannon_evenness(as.numeric(table(x)))
}

#' Fraction accounted for by the most frequent category
#'
#' @param x Character vector.
#' @return Numeric in (0, 1], or `NA` if empty.
top_fraction <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_real_)
  max(table(x)) / length(x)
}

#' Publicness of an epitope's repertoire
#'
#' A CDR3 seen in several independent donors is "public": easy to generate and
#' commonly selected. Publicness is computed only over records that carry a
#' donor ID, and is `NA` when fewer than `min_donors` donors are available -
#' with one donor, nothing can be shared, and reporting 0 would be misleading.
#'
#' @param cdr3 Character vector of CDR3 sequences.
#' @param donor Character vector of donor keys, same length as `cdr3`.
#' @param min_donors Minimum donors required to report publicness.
#' @return A list of publicness statistics.
publicness <- function(cdr3, donor, min_donors = 2L) {
  ok <- !is.na(cdr3) & !is.na(donor)
  cdr3 <- cdr3[ok]; donor <- donor[ok]

  n_donors <- length(unique(donor))
  if (n_donors < min_donors) {
    return(list(
      n_donors_used       = n_donors,
      public_fraction     = NA_real_,
      mean_donors_per_cdr3 = NA_real_,
      max_donors_per_cdr3  = NA_integer_
    ))
  }

  dt <- data.table::data.table(cdr3 = cdr3, donor = donor)
  per_cdr3 <- dt[, .(n_donors = data.table::uniqueN(donor)), by = cdr3]

  list(
    n_donors_used        = n_donors,
    public_fraction      = mean(per_cdr3$n_donors >= 2),
    mean_donors_per_cdr3 = mean(per_cdr3$n_donors),
    max_donors_per_cdr3  = max(per_cdr3$n_donors)
  )
}

#' Compute the per-epitope feature table
#'
#' @param dt A table from `scripts/02_parse_vdjdb.R` (one row per TCR).
#' @param rarefy_depth Depth for rarefied richness. Epitopes with fewer records
#'   get `NA`, which is deliberate: an un-comparable value is better than a
#'   silently biased one.
#' @param min_tcr Minimum records for an epitope to be included at all.
#' @return A `data.table`, one row per epitope.
epitope_features <- function(dt, rarefy_depth = 30L, min_tcr = 10L) {
  stopifnot(all(c("epitope", "cdr3") %in% names(dt)))

  has_donor <- "donor_key" %in% names(dt)
  if (!has_donor) {
    log_msg("no donor_key column: publicness features will be NA", level = "warn")
  }

  feats <- dt[, {
    cdr3_counts <- as.numeric(table(cdr3))
    pub <- if (has_donor) publicness(cdr3, donor_key) else
      list(n_donors_used = NA_integer_, public_fraction = NA_real_,
           mean_donors_per_cdr3 = NA_real_, max_donors_per_cdr3 = NA_integer_)
    mhcs <- sort(unique(mhc_a[!is.na(mhc_a)]))

    .(
      ## --- scale -------------------------------------------------------------
      n_tcr          = .N,
      n_cdr3         = data.table::uniqueN(cdr3),
      n_studies      = data.table::uniqueN(reference_id),
      n_donors       = pub$n_donors_used,

      ## --- publicness (the repertoire-availability proxy) ---------------------
      public_fraction      = pub$public_fraction,
      mean_donors_per_cdr3 = pub$mean_donors_per_cdr3,
      max_donors_per_cdr3  = pub$max_donors_per_cdr3,

      ## --- diversity, raw and rarefied ----------------------------------------
      shannon        = shannon_entropy(cdr3_counts),
      evenness       = shannon_evenness(cdr3_counts),
      clonality      = clonality(cdr3_counts),
      inv_simpson    = inverse_simpson(cdr3_counts),
      chao1          = chao1(cdr3_counts),
      richness_rare  = rarefy_richness(cdr3_counts, rarefy_depth),

      ## --- sequence structure --------------------------------------------------
      cdr3_len_mean  = mean(nchar(cdr3)),
      cdr3_len_sd    = stats::sd(nchar(cdr3)),

      ## --- germline usage ------------------------------------------------------
      v_entropy      = usage_entropy(v_segm),
      j_entropy      = usage_entropy(j_segm),
      v_top_fraction = top_fraction(v_segm),

      ## --- provenance ----------------------------------------------------------
      pct_direct     = 100 * mean(direct_binding, na.rm = TRUE),
      mhc_class      = mhc_class[1],
      mhc_a          = if (length(mhcs)) mhcs[1] else NA_character_,
      n_mhc_a        = length(mhcs),
      protein        = protein_label[1],
      organism       = organism[1],
      structural     = structural[1]
    )
  }, by = epitope]

  out <- feats[n_tcr >= min_tcr]
  log_msg("features: %d epitopes retained of %d (min_tcr = %d)",
          nrow(out), nrow(feats), min_tcr)
  log_msg("  with rarefied richness at depth %d: %d",
          rarefy_depth, sum(!is.na(out$richness_rare)))
  log_msg("  with publicness (>= 2 donors): %d",
          sum(!is.na(out$public_fraction)))
  out[order(-n_donors, -n_tcr)]
}

#' Restrict to epitopes where the metrics are actually interpretable
#'
#' Publicness needs several donors; rarefied richness needs enough records.
#' Rather than quietly comparing incomparable rows, this returns the analysable
#' subset and logs what was dropped and why.
#'
#' @param feats Output of `epitope_features()`.
#' @param min_donors Minimum donors.
#' @param min_tcr Minimum records.
#' @param require_rarefied Drop epitopes whose rarefied richness is `NA`.
#' @return A filtered `data.table`.
analysable_epitopes <- function(feats, min_donors = 3L, min_tcr = 20L,
                                require_rarefied = TRUE) {
  n0 <- nrow(feats)
  out <- feats[!is.na(n_donors) & n_donors >= min_donors & n_tcr >= min_tcr]
  if (require_rarefied) out <- out[!is.na(richness_rare)]

  log_msg("analysable set: %d of %d epitopes (>= %d donors, >= %d records)",
          nrow(out), n0, min_donors, min_tcr)
  if (nrow(out)) {
    log_msg("  by protein: %s",
            paste(sprintf("%s=%d", names(table(out$protein)),
                          as.integer(table(out$protein))), collapse = ", "))
  }
  out[]
}
