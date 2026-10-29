## R/cdr3_network.R ----------------------------------------------------------
## CDR3 similarity networks, and publicness measured over clusters.
##
## WHY THIS EXISTS
##
## Exact-match publicness is near-zero in VDJdb: across the analysable
## SARS-CoV-2 structural epitopes, the mean number of donors per unique CDR3 is
## about 1.0-1.17, and most epitopes have a public fraction of 0.
##
## That is not evidence that these repertoires are private. T cells responding
## to the same epitope in different people converge on shared *motifs* rather
## than on identical sequences, and each donor contributes only a handful of
## records, so exact matches are rare by construction.
##
## Measuring sharing over clusters of similar CDR3s rather than over identical
## strings is therefore not a refinement - it is what makes publicness
## measurable at all in this dataset. Both versions are reported side by side
## so the difference is visible rather than assumed.
## -----------------------------------------------------------------------------

#' Cluster CDR3 sequences by edit distance
#'
#' Builds a graph linking sequences within `max_dist` edits and takes connected
#' components as clusters. Single-linkage is deliberate: it is the standard
#' construction for TCR similarity networks and is permissive enough to join a
#' motif family, but it can chain, so `largest_cluster_fraction` is reported as
#' a diagnostic.
#'
#' @param cdr3 Character vector of CDR3 sequences (duplicates allowed).
#' @param max_dist Maximum edit distance for an edge. 1 is conservative; 2 is
#'   common in the literature and far more permissive.
#' @param method Passed to `stringdist` ("lv" = Levenshtein, "hamming" =
#'   substitutions only, which forces equal lengths).
#' @param max_unique Guard against quadratic blow-up. Above this many unique
#'   sequences the function returns `NULL` rather than attempting the matrix.
#' @return A `data.table` with `cdr3` and `cluster`, one row per unique
#'   sequence, or `NULL` if the input is too large or too small.
cdr3_clusters <- function(cdr3, max_dist = 1L, method = "lv",
                          max_unique = 4000L) {
  if (!requireNamespace("stringdist", quietly = TRUE) ||
      !requireNamespace("igraph", quietly = TRUE)) {
    stop("cdr3_clusters() needs the stringdist and igraph packages.", call. = FALSE)
  }

  u <- unique(cdr3[!is.na(cdr3) & nzchar(cdr3)])
  n <- length(u)
  if (n < 2L) return(NULL)
  if (n > max_unique) {
    log_msg("skipping network: %d unique CDR3s exceeds max_unique = %d",
            n, max_unique, level = "warn")
    return(NULL)
  }

  dm <- stringdist::stringdistmatrix(u, u, method = method)
  adj <- dm <= max_dist
  diag(adj) <- FALSE

  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  igraph::V(g)$name <- u
  comp <- igraph::components(g)

  data.table::data.table(cdr3 = u, cluster = as.integer(comp$membership))
}

#' Publicness measured over CDR3 clusters rather than exact sequences
#'
#' @param cdr3 Character vector of CDR3 sequences.
#' @param donor Character vector of donor keys, same length.
#' @param max_dist Edit-distance threshold for clustering.
#' @param min_donors Minimum donors required to report publicness.
#' @return A named list of cluster-level statistics.
cluster_publicness <- function(cdr3, donor, max_dist = 1L, min_donors = 2L) {
  na_result <- list(
    n_clusters                 = NA_integer_,
    cluster_public_fraction    = NA_real_,
    mean_donors_per_cluster    = NA_real_,
    max_donors_per_cluster     = NA_integer_,
    largest_cluster_fraction   = NA_real_,
    clustered_fraction         = NA_real_
  )

  ok <- !is.na(cdr3) & !is.na(donor)
  cdr3 <- cdr3[ok]; donor <- donor[ok]
  if (length(unique(donor)) < min_donors) return(na_result)

  cl <- cdr3_clusters(cdr3, max_dist = max_dist)
  if (is.null(cl)) return(na_result)

  dt <- merge(data.table::data.table(cdr3 = cdr3, donor = donor),
              cl, by = "cdr3", all.x = TRUE)
  dt <- dt[!is.na(cluster)]
  if (!nrow(dt)) return(na_result)

  per_cluster <- dt[, .(n_donors = data.table::uniqueN(donor),
                        n_cdr3   = data.table::uniqueN(cdr3)), by = cluster]

  list(
    n_clusters               = nrow(per_cluster),
    cluster_public_fraction  = mean(per_cluster$n_donors >= 2),
    mean_donors_per_cluster  = mean(per_cluster$n_donors),
    max_donors_per_cluster   = max(per_cluster$n_donors),
    ## Chaining diagnostic: if one component swallows most sequences, the
    ## clustering is too permissive and the publicness figure is inflated.
    largest_cluster_fraction = max(per_cluster$n_cdr3) / sum(per_cluster$n_cdr3),
    ## Fraction of sequences with at least one neighbour - a direct measure of
    ## how structured the epitope's repertoire is.
    clustered_fraction       = mean(per_cluster$n_cdr3[match(dt$cluster, per_cluster$cluster)] > 1)
  )
}

#' Add cluster-based features to an existing per-epitope table
#'
#' @param feats Output of `epitope_features()`.
#' @param dt The TCR-level table those features were built from.
#' @param max_dist Edit-distance threshold.
#' @return `feats` with cluster columns joined on `epitope`.
add_cluster_features <- function(feats, dt, max_dist = 1L) {
  eps <- feats$epitope
  log_msg("building CDR3 networks for %d epitopes (max_dist = %d)",
          length(eps), max_dist)

  rows <- lapply(eps, function(ep) {
    sub <- dt[epitope == ep]
    res <- cluster_publicness(sub$cdr3, sub$donor_key, max_dist = max_dist)
    data.table::as.data.table(c(list(epitope = ep), res))
  })

  clust <- data.table::rbindlist(rows, fill = TRUE)
  out <- merge(feats, clust, by = "epitope", all.x = TRUE, sort = FALSE)

  n_ok <- sum(!is.na(out$cluster_public_fraction))
  log_msg("  cluster publicness computed for %d of %d epitopes", n_ok, nrow(out))
  out[order(-n_donors, -n_tcr)]
}
