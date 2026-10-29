## R/diversity.R -------------------------------------------------------------
## Diversity, evenness and richness estimators.
##
## Implemented directly rather than via immunarch: they are short, they make
## the definitions explicit in the repository, and they avoid compiling a CRAN
## package against a conda R.
##
## IMPORTANT: every index here is computed on *counts of unique CDR3s within an
## epitope's record set*. These are not clone sizes. VDJdb records reflect
## study effort, not cell numbers, so "diversity" here means diversity of the
## sequences reported against an epitope - a property of the repertoire that
## recognises it, not of any individual's clonal expansion.
## -----------------------------------------------------------------------------

#' Shannon entropy
#'
#' @param counts Numeric vector of counts (zeros and NAs dropped).
#' @param base Logarithm base; natural log by default.
#' @return Entropy in nats (or the chosen base). `NA` for an empty vector.
shannon_entropy <- function(counts, base = exp(1)) {
  counts <- counts[!is.na(counts) & counts > 0]
  if (!length(counts)) return(NA_real_)
  p <- counts / sum(counts)
  -sum(p * log(p, base = base))
}

#' Normalised Shannon evenness, in [0, 1]
#'
#' Entropy divided by its maximum for the observed richness. Defined as 1 for a
#' single type, since a one-type sample is perfectly even.
#'
#' @param counts Numeric vector of counts.
#' @return Evenness in [0, 1].
shannon_evenness <- function(counts) {
  counts <- counts[!is.na(counts) & counts > 0]
  n <- length(counts)
  if (!n) return(NA_real_)
  if (n == 1L) return(1)
  shannon_entropy(counts) / log(n)
}

#' Clonality: 1 - evenness
#'
#' High clonality means a few sequences dominate the set.
#'
#' @param counts Numeric vector of counts.
#' @return Clonality in [0, 1].
clonality <- function(counts) {
  e <- shannon_evenness(counts)
  if (is.na(e)) return(NA_real_)
  1 - e
}

#' Simpson's index and its inverse
#'
#' `simpson_index()` is the probability that two draws are the same type.
#' `inverse_simpson()` is its reciprocal, an effective number of types.
#'
#' @param counts Numeric vector of counts.
#' @return Numeric scalar.
simpson_index <- function(counts) {
  counts <- counts[!is.na(counts) & counts > 0]
  if (!length(counts)) return(NA_real_)
  p <- counts / sum(counts)
  sum(p^2)
}

#' @rdname simpson_index
inverse_simpson <- function(counts) {
  s <- simpson_index(counts)
  if (is.na(s) || s == 0) return(NA_real_)
  1 / s
}

#' Chao1 richness estimator (bias-corrected)
#'
#' Estimates how many types exist, including those not seen, from the number of
#' singletons and doubletons. With no singletons it returns observed richness.
#'
#' Treat Chao1 cautiously here: it assumes random sampling from a fixed
#' community, and TCR records accumulate through study effort instead. It is
#' reported as a descriptive lower bound, not a true repertoire size.
#'
#' @param counts Numeric vector of counts.
#' @return Estimated richness.
chao1 <- function(counts) {
  counts <- counts[!is.na(counts) & counts > 0]
  if (!length(counts)) return(NA_real_)
  s_obs <- length(counts)
  f1 <- sum(counts == 1)
  f2 <- sum(counts == 2)
  s_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))
}

#' Hill number of order q
#'
#' Unifies the indices above: q = 0 is richness, q -> 1 is exp(Shannon),
#' q = 2 is inverse Simpson. Larger q weights common types more heavily.
#'
#' @param counts Numeric vector of counts.
#' @param q Order (>= 0).
#' @return Effective number of types.
hill_number <- function(counts, q = 1) {
  counts <- counts[!is.na(counts) & counts > 0]
  if (!length(counts)) return(NA_real_)
  p <- counts / sum(counts)
  if (abs(q - 1) < 1e-8) return(exp(shannon_entropy(counts)))
  if (q == 0) return(length(p))
  sum(p^q)^(1 / (1 - q))
}

#' Expected richness at a reduced sampling depth (rarefaction)
#'
#' Uses Hurlbert's analytic expectation rather than repeated subsampling, so
#' the result is exact and deterministic.
#'
#' Rarefaction is essential before comparing epitopes: richness grows with the
#' number of records, and record counts here vary by three orders of magnitude.
#'
#' @param counts Numeric vector of counts.
#' @param depth Number of records to rarefy to.
#' @return Expected number of distinct types at `depth`, or `NA` if the sample
#'   is smaller than `depth`.
rarefy_richness <- function(counts, depth) {
  counts <- counts[!is.na(counts) & counts > 0]
  n <- sum(counts)
  if (!length(counts) || depth > n) return(NA_real_)
  if (depth == n) return(as.numeric(length(counts)))
  ## E[S] = sum_i (1 - C(n - c_i, m) / C(n, m)), computed in log space.
  log_denom <- lchoose(n, depth)
  terms <- vapply(counts, function(ci) {
    if (n - ci < depth) return(1)
    1 - exp(lchoose(n - ci, depth) - log_denom)
  }, numeric(1))
  sum(terms)
}

#' Repertoire overlap between two sets of sequences
#'
#' @param a,b Character vectors of sequences (duplicates allowed).
#' @param method "jaccard" (presence/absence) or "morisita" (abundance-weighted).
#' @return Overlap in [0, 1], or `NA` if either set is empty.
repertoire_overlap <- function(a, b, method = c("jaccard", "morisita")) {
  method <- match.arg(method)
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  if (!length(a) || !length(b)) return(NA_real_)

  if (method == "jaccard") {
    ua <- unique(a); ub <- unique(b)
    return(length(intersect(ua, ub)) / length(union(ua, ub)))
  }

  ta <- table(a); tb <- table(b)
  shared <- intersect(names(ta), names(tb))
  if (!length(shared)) return(0)
  na <- sum(ta); nb <- sum(tb)
  xa <- as.numeric(ta[shared]); xb <- as.numeric(tb[shared])
  num <- 2 * sum(xa * xb) / (na * nb)
  den <- simpson_index(as.numeric(ta)) + simpson_index(as.numeric(tb))
  if (den == 0) return(NA_real_)
  num / den
}
