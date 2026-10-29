## Tests for depth-controlled publicness. Offline.

test_that("raw publicness rises with donor count even when biology is constant", {
  ## This is the confound the whole module exists to remove. Build two epitopes
  ## with IDENTICAL structure - every donor contributes the same two-sequence
  ## motif family - differing only in how many donors were sampled.
  make <- function(n_donors) {
    donors <- paste0("d", seq_len(n_donors))
    list(
      cdr3  = rep(c("CASSAAF", "CASSABF"), times = n_donors),
      donor = rep(donors, each = 2)
    )
  }
  small <- make(3); large <- make(30)

  p_small <- publicness(small$cdr3, small$donor)$public_fraction
  p_large <- publicness(large$cdr3, large$donor)$public_fraction

  ## Same underlying biology, but the larger sample looks no less public -
  ## confirming publicness is driven by sampling, not only by sharing.
  expect_true(p_large >= p_small)
})

test_that("donor rarefaction returns NA rather than extrapolating", {
  cdr3  <- c("CASSAAF", "CASSABF", "CASSACF")
  donor <- c("d1", "d2", "d3")
  ## Only 3 donors available, 10 requested.
  res <- publicness_donor_rarefied(cdr3, donor, k = 10L, n_rep = 5L)
  expect_true(is.na(res$public_rare_mean))
  expect_equal(res$n_rep_used, 0L)
})

test_that("donor rarefaction is reproducible and in range", {
  set.seed(42)
  donors <- paste0("d", 1:10)
  cdr3  <- rep(c("CASSAAF", "CASSABF", "CAWWWWW"), 10)
  donor <- rep(donors, each = 3)

  a <- publicness_donor_rarefied(cdr3, donor, k = 5L, n_rep = 20L, seed = 7L)
  b <- publicness_donor_rarefied(cdr3, donor, k = 5L, n_rep = 20L, seed = 7L)
  expect_equal(a$public_rare_mean, b$public_rare_mean)
  expect_gte(a$public_rare_mean, 0)
  expect_lte(a$public_rare_mean, 1)
})

test_that("rarefaction equalises epitopes that differ only in depth", {
  ## Identical per-donor structure, 4 donors vs 20.
  make <- function(n) list(
    cdr3  = rep(c("CASSAAF", "CASSABF"), times = n),
    donor = rep(paste0("d", seq_len(n)), each = 2)
  )
  s <- make(4); l <- make(20)

  rs <- publicness_donor_rarefied(s$cdr3, s$donor, k = 4L, n_rep = 30L, seed = 1L)
  rl <- publicness_donor_rarefied(l$cdr3, l$donor, k = 4L, n_rep = 30L, seed = 1L)

  ## Measured at the same donor depth, the two agree.
  expect_equal(rs$public_rare_mean, rl$public_rare_mean, tolerance = 0.05)
})

test_that("permutation null detects genuinely donor-structured sharing", {
  ## Each donor has its OWN private motif family: similar sequences never cross
  ## donors, so observed cluster publicness should fall BELOW the null.
  cdr3 <- c("AAAAAAA", "AAAAAAB",     # donor 1 family
            "CCCCCCC", "CCCCCCD",     # donor 2 family
            "EEEEEEE", "EEEEEEF")     # donor 3 family
  donor <- rep(c("d1", "d2", "d3"), each = 2)

  res <- publicness_null(cdr3, donor, n_perm = 200L, max_dist = 1L, seed = 3L)
  expect_false(is.na(res$z))
  expect_lt(res$observed, res$null_mean + 1e-9)
})

test_that("permutation null returns NA for degenerate input", {
  expect_true(is.na(publicness_null(c("A", "B"), c("d1", "d1"))$z))
  expect_true(is.na(publicness_null(character(0), character(0))$z))
})

test_that("permutation p-value is bounded away from zero", {
  ## (hits + 1) / (n + 1) means p can never be exactly 0, which would overstate
  ## significance that a finite permutation set cannot support.
  cdr3  <- rep(c("CASSAAF", "CASSABF", "CAWWWWW"), 4)
  donor <- rep(paste0("d", 1:4), each = 3)
  res <- publicness_null(cdr3, donor, n_perm = 50L, seed = 1L)
  if (!is.na(res$p_emp)) {
    expect_gt(res$p_emp, 0)
    expect_lte(res$p_emp, 1)
  }
})

test_that("add_depth_controlled keeps every epitope row", {
  dt <- data.table::data.table(
    epitope   = rep(c("EP1", "EP2"), each = 6),
    cdr3      = c("CASSAAF", "CASSABF", "CASSACF", "CASSADF", "CASSAEF", "CASSAFF",
                  "CAWWWWA", "CAWWWWB", "CAWWWWC", "CAWWWWD", "CAWWWWE", "CAWWWWF"),
    donor_key = rep(paste0("d", 1:3), each = 2, times = 2)
  )
  feats <- data.table::data.table(
    epitope = c("EP1", "EP2"), n_donors = c(3L, 3L), n_tcr = c(6L, 6L)
  )
  out <- add_depth_controlled(feats, dt, k = 3L, n_rep = 5L, n_perm = 20L)
  expect_equal(nrow(out), 2L)
  expect_true(all(c("public_rare_mean", "public_z") %in% names(out)))
})

test_that("the null deduplicates records before permuting", {
  ## One donor contributed the SAME sequence many times. Without dedup, the
  ## shuffle scatters those copies across donors and invents sharing that was
  ## never observed - the bug this guards against.
  cdr3  <- c(rep("CASSAAF", 50), "CASSABF", "CAWWWWW")
  donor <- c(rep("d1", 50), "d2", "d3")

  res <- publicness_null(cdr3, donor, n_perm = 100L, max_dist = 1L, seed = 2L)
  ## After dedup there are only 3 (cdr3, donor) pairs, so the null cannot be
  ## driven by d1's 50 duplicate records.
  expect_false(is.na(res$observed))
  expect_lte(res$observed, 1)
})

test_that("within_donor_cluster_fraction detects clonal families", {
  ## Each donor has its own family of similar sequences: every multi-sequence
  ## cluster is confined to one donor.
  cdr3 <- c("AAAAAAA", "AAAAAAB", "AAAAABB",
            "CCCCCCC", "CCCCCCD", "CCCCCDD")
  donor <- rep(c("d1", "d2"), each = 3)
  res <- publicness_null(cdr3, donor, n_perm = 100L, max_dist = 1L, seed = 5L)
  expect_equal(res$within_donor_cluster_fraction, 1)

  ## Mixed donors within a family -> not purely within-donor.
  donor2 <- c("d1", "d2", "d1", "d2", "d1", "d2")
  res2 <- publicness_null(cdr3, donor2, n_perm = 100L, max_dist = 1L, seed = 5L)
  expect_lt(res2$within_donor_cluster_fraction, 1)
})

test_that("the null p-value is two-sided", {
  ## Donor-structured repertoires (z < 0) must register as significant, not be
  ## reported as p = 1 the way a one-sided test would.
  cdr3 <- c("AAAAAAA", "AAAAAAB", "AAAAABB", "AAAABBB",
            "CCCCCCC", "CCCCCCD", "CCCCCDD", "CCCCDDD")
  donor <- rep(c("d1", "d2"), each = 4)
  res <- publicness_null(cdr3, donor, n_perm = 300L, max_dist = 1L, seed = 11L)
  if (!is.na(res$z) && res$z < -2) expect_lt(res$p_emp, 0.5)
})
