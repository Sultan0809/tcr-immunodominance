## Tests for diversity estimators and per-epitope features. Fully offline.

test_that("shannon_entropy matches known values", {
  ## Four equally abundant types: H = log(4).
  expect_equal(shannon_entropy(c(1, 1, 1, 1)), log(4))
  ## A single type carries no information.
  expect_equal(shannon_entropy(c(5)), 0)
  ## Base-2 entropy of a fair 2-way split is 1 bit.
  expect_equal(shannon_entropy(c(1, 1), base = 2), 1)
  ## Zeros and NAs must not contribute.
  expect_equal(shannon_entropy(c(1, 1, 0, NA)), log(2))
  expect_true(is.na(shannon_entropy(numeric(0))))
})

test_that("evenness and clonality are complementary and bounded", {
  even <- c(10, 10, 10, 10)
  skew <- c(97, 1, 1, 1)

  expect_equal(shannon_evenness(even), 1)
  expect_equal(clonality(even), 0)
  expect_lt(shannon_evenness(skew), 0.5)
  expect_gt(clonality(skew), 0.5)

  ## A single type is defined as perfectly even.
  expect_equal(shannon_evenness(c(7)), 1)
  expect_equal(clonality(c(7)), 0)

  for (x in list(even, skew, c(1, 2, 3))) {
    expect_gte(shannon_evenness(x), 0)
    expect_lte(shannon_evenness(x), 1)
  }
})

test_that("Simpson behaves as an effective number of types", {
  expect_equal(simpson_index(c(1, 1)), 0.5)
  expect_equal(inverse_simpson(c(1, 1, 1, 1)), 4)
  expect_equal(inverse_simpson(c(1)), 1)
})

test_that("chao1 adds unseen types only when singletons exist", {
  ## No singletons: estimate equals observed richness.
  expect_equal(chao1(c(2, 2, 3)), 3)
  ## Singletons present: estimate exceeds observed.
  expect_gt(chao1(c(1, 1, 1, 5)), 4)
  expect_true(is.na(chao1(numeric(0))))
})

test_that("hill numbers reduce to the standard indices", {
  x <- c(5, 3, 2, 1)
  expect_equal(hill_number(x, 0), 4)
  expect_equal(hill_number(x, 1), exp(shannon_entropy(x)), tolerance = 1e-8)
  expect_equal(hill_number(x, 2), inverse_simpson(x), tolerance = 1e-8)
  ## Hill numbers are non-increasing in q.
  expect_gte(hill_number(x, 0), hill_number(x, 1))
  expect_gte(hill_number(x, 1), hill_number(x, 2))
})

test_that("rarefy_richness interpolates and refuses extrapolation", {
  x <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)   # 10 singletons, n = 10
  ## At full depth, expected richness is observed richness.
  expect_equal(rarefy_richness(x, 10), 10)
  ## At half depth, exactly half the types are expected to appear.
  expect_equal(rarefy_richness(x, 5), 5, tolerance = 1e-8)
  ## Extrapolation beyond the sample is refused rather than guessed.
  expect_true(is.na(rarefy_richness(x, 20)))

  ## Rarefaction must be monotonic in depth.
  y <- c(50, 20, 10, 5, 5, 3, 2, 1, 1, 1)
  depths <- c(5, 10, 25, 50)
  vals <- vapply(depths, function(d) rarefy_richness(y, d), numeric(1))
  expect_true(all(diff(vals) >= 0))
})

test_that("rarefaction removes the sample-size confound it is there to remove", {
  ## Two epitopes with identical underlying diversity but 10x different depth.
  set.seed(1)
  pool <- paste0("CDR3_", 1:40)
  small <- sample(pool, 40, replace = TRUE)
  large <- sample(pool, 400, replace = TRUE)

  ## Raw richness differs substantially...
  expect_gt(length(unique(large)), length(unique(small)))
  ## ...but rarefied to a common depth the two agree closely.
  r_small <- rarefy_richness(as.numeric(table(small)), 30)
  r_large <- rarefy_richness(as.numeric(table(large)), 30)
  expect_lt(abs(r_small - r_large), 5)
})

test_that("repertoire_overlap handles identical, disjoint and empty sets", {
  a <- c("X", "Y", "Z")
  expect_equal(repertoire_overlap(a, a, "jaccard"), 1)
  expect_equal(repertoire_overlap(a, c("P", "Q"), "jaccard"), 0)
  expect_equal(repertoire_overlap(a, c("Y", "Z", "W"), "jaccard"), 2 / 4)
  expect_true(is.na(repertoire_overlap(a, character(0))))

  expect_equal(repertoire_overlap(a, a, "morisita"), 1, tolerance = 1e-8)
  expect_equal(repertoire_overlap(a, c("P"), "morisita"), 0)
})

test_that("usage_entropy and top_fraction describe germline skew", {
  expect_equal(usage_entropy(c("V1", "V2", "V3", "V4")), 1)
  expect_lt(usage_entropy(c("V1", "V1", "V1", "V2")), 1)
  expect_equal(top_fraction(c("V1", "V1", "V1", "V2")), 0.75)
  expect_true(is.na(usage_entropy(c(NA, ""))))
})

test_that("publicness is NA with one donor, not zero", {
  ## With a single donor nothing CAN be shared; reporting 0 would read as
  ## "this epitope's repertoire is private", which the data cannot support.
  res <- publicness(c("A", "B", "C"), c("d1", "d1", "d1"))
  expect_equal(res$n_donors_used, 1L)
  expect_true(is.na(res$public_fraction))
})

test_that("publicness counts cross-donor sharing correctly", {
  cdr3  <- c("A", "A", "B", "C")
  donor <- c("d1", "d2", "d1", "d2")
  res <- publicness(cdr3, donor)

  expect_equal(res$n_donors_used, 2L)
  ## A is in both donors; B and C in one each -> 1 of 3 unique CDR3s is public.
  expect_equal(res$public_fraction, 1 / 3)
  expect_equal(res$max_donors_per_cdr3, 2L)
})

test_that("publicness ignores records without a donor ID", {
  res <- publicness(c("A", "A", "B"), c("d1", NA, "d2"))
  expect_equal(res$n_donors_used, 2L)
  ## A appears once with a donor, so nothing is shared across donors.
  expect_equal(res$public_fraction, 0)
})

test_that("epitope_features returns one row per epitope and respects min_tcr", {
  dt <- data.table::data.table(
    epitope        = rep(c("EP1", "EP2"), times = c(12, 3)),
    cdr3           = c(paste0("CASS", 1:12), paste0("CASW", 1:3)),
    v_segm         = rep("TRBV7-9", 15),
    j_segm         = rep("TRBJ2-7", 15),
    mhc_a          = "HLA-A*02:01",
    mhc_class      = "MHCI",
    donor_key      = rep(c("s|d1", "s|d2", "s|d3"), 5),
    reference_id   = "PMID:1",
    direct_binding = TRUE,
    protein_label  = "S",
    organism       = "SARS-CoV-2",
    structural     = TRUE
  )

  f <- epitope_features(dt, rarefy_depth = 5L, min_tcr = 10L)
  expect_equal(nrow(f), 1L)          # EP2 has only 3 records
  expect_equal(f$epitope, "EP1")
  expect_equal(f$n_tcr, 12L)
  expect_equal(f$n_donors, 3L)

  ## Depth beyond the sample yields NA rather than a biased number.
  f2 <- epitope_features(dt, rarefy_depth = 500L, min_tcr = 1L)
  expect_true(all(is.na(f2$richness_rare)))
})

test_that("analysable_epitopes drops uninterpretable rows", {
  feats <- data.table::data.table(
    epitope       = c("A", "B", "C"),
    n_donors      = c(10L, 1L, 5L),
    n_tcr         = c(100L, 50L, 5L),
    richness_rare = c(20, NA, 4),
    protein       = c("S", "S", "N")
  )
  out <- analysable_epitopes(feats, min_donors = 3L, min_tcr = 20L)
  expect_equal(out$epitope, "A")     # B too few donors, C too few records
})
