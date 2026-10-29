## Tests for CDR3 similarity networks and cluster publicness. Offline.

test_that("cdr3_clusters groups sequences within the distance threshold", {
  ## Two families: CASSAA.. variants and a distant CAWXYZ.
  seqs <- c("CASSAAF", "CASSABF", "CASSACF", "CAWWWWW")
  cl <- cdr3_clusters(seqs, max_dist = 1L)

  expect_equal(nrow(cl), 4L)
  ## The three one-edit-apart sequences share a cluster.
  fam <- cl[cdr3 %in% c("CASSAAF", "CASSABF", "CASSACF"), cluster]
  expect_equal(length(unique(fam)), 1L)
  ## The distant sequence is on its own.
  expect_false(cl[cdr3 == "CAWWWWW", cluster] %in% fam)
})

test_that("a stricter threshold yields more clusters", {
  seqs <- c("CASSAAF", "CASSABF", "CASSBBF", "CASSCCF")
  n1 <- data.table::uniqueN(cdr3_clusters(seqs, max_dist = 1L)$cluster)
  n2 <- data.table::uniqueN(cdr3_clusters(seqs, max_dist = 2L)$cluster)
  expect_gte(n1, n2)
})

test_that("cdr3_clusters refuses degenerate and oversized input", {
  expect_null(cdr3_clusters(character(0)))
  expect_null(cdr3_clusters("CASSAAF"))            # single sequence
  expect_null(cdr3_clusters(c(NA, NA)))
  ## Oversized input returns NULL rather than attempting a huge matrix.
  big <- paste0("CASS", sprintf("%05d", 1:50))
  expect_null(cdr3_clusters(big, max_unique = 10L))
})

test_that("cluster publicness finds sharing that exact matching misses", {
  ## Two donors with SIMILAR but never IDENTICAL CDR3s - the situation that
  ## makes exact-match publicness read as zero in VDJdb.
  cdr3  <- c("CASSAAF", "CASSABF", "CASSACF", "CASSADF")
  donor <- c("d1",      "d1",      "d2",      "d2")

  exact <- publicness(cdr3, donor)
  expect_equal(exact$public_fraction, 0)      # no sequence is shared outright

  clust <- cluster_publicness(cdr3, donor, max_dist = 1L)
  expect_equal(clust$n_clusters, 1L)
  expect_equal(clust$cluster_public_fraction, 1)   # the cluster spans both donors
  expect_equal(clust$max_donors_per_cluster, 2L)
})

test_that("cluster publicness stays NA with a single donor", {
  res <- cluster_publicness(c("CASSAAF", "CASSABF"), c("d1", "d1"))
  expect_true(is.na(res$cluster_public_fraction))
  expect_true(is.na(res$n_clusters))
})

test_that("genuinely unrelated repertoires are not called public", {
  cdr3  <- c("CASSAAF", "CAWWWWW")
  donor <- c("d1", "d2")
  res <- cluster_publicness(cdr3, donor, max_dist = 1L)
  expect_equal(res$n_clusters, 2L)
  expect_equal(res$cluster_public_fraction, 0)
})

test_that("largest_cluster_fraction flags single-linkage chaining", {
  ## A chain where each sequence is one edit from the next: single linkage
  ## merges them all, which is exactly the failure mode to detect.
  chain <- c("AAAAAAA", "AAAAAAB", "AAAAABB", "AAAABBB", "AAABBBB")
  res <- cluster_publicness(chain, rep(c("d1", "d2"), length.out = 5),
                            max_dist = 1L)
  expect_equal(res$n_clusters, 1L)
  expect_equal(res$largest_cluster_fraction, 1)
})

test_that("add_cluster_features joins without dropping epitopes", {
  dt <- data.table::data.table(
    epitope   = rep(c("EP1", "EP2"), each = 4),
    cdr3      = c("CASSAAF", "CASSABF", "CASSACF", "CASSADF",
                  "CAWWWWA", "CAWWWWB", "CAWWWWC", "CAWWWWD"),
    donor_key = rep(c("d1", "d1", "d2", "d2"), 2)
  )
  feats <- data.table::data.table(
    epitope = c("EP1", "EP2"), n_donors = c(2L, 2L), n_tcr = c(4L, 4L)
  )
  out <- add_cluster_features(feats, dt, max_dist = 1L)
  expect_equal(nrow(out), 2L)
  expect_setequal(out$epitope, c("EP1", "EP2"))
  expect_true(all(out$cluster_public_fraction == 1))
})
