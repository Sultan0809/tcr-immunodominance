## Tests for predicted-epitope matching. Offline.

test_that("lcs_length finds the longest shared run", {
  expect_equal(lcs_length("ABCDEF", "XXCDEXX")$len, 3L)
  expect_equal(lcs_length("ABCDEF", "XXCDEXX")$seq, "CDE")
  expect_equal(lcs_length("ABC", "ABC")$len, 3L)
  expect_equal(lcs_length("ABC", "XYZ")$len, 0L)
  expect_equal(lcs_length("", "ABC")$len, 0L)
  expect_equal(lcs_length(NA_character_, "ABC")$len, 0L)
})

test_that("classify_match distinguishes containment from partial overlap", {
  expect_equal(classify_match("YLQPRTFLL", "YLQPRTFLL")$type, "exact")

  ## A class II 15-mer containing the database's 14-mer.
  expect_equal(classify_match("TFEYVSQPFLMDLEG", "TFEYVSQPFLMDLE")$type,
               "db within predicted")

  ## A 9-mer inside a longer database entry.
  expect_equal(classify_match("PITVYYAVL", "PPSIPITVYYAVLERA")$type,
               "predicted within db")

  ## Offset overlap: neither contains the other, but they share a core.
  m <- classify_match("SLPITVYYA", "PPSIPITVYYAVLERA", min_core = 7L)
  expect_equal(m$type, "shared core")
  expect_equal(m$core, "PITVYYA")

  ## Unrelated sequences must not match.
  expect_true(is.na(classify_match("YLQPRTFLL", "KKKKKKKKK")$type))
})

test_that("min_core controls match stringency", {
  ## A 5-residue overlap passes at min_core 5 and fails at 7.
  expect_equal(classify_match("AAAAACDEFG", "ZZZZZCDEFG", min_core = 5L)$type,
               "shared core")
  a <- classify_match("ABCDEWWWW", "ZZZZZABCDE", min_core = 7L)
  expect_true(is.na(a$type))
})

test_that("unmatched predicted epitopes are retained, not dropped", {
  predicted <- data.table::data.table(
    epitope = c("YLQPRTFLL", "KKKKKKKKK"),
    organism = "SARS-CoV-2", mhc_class = "I"
  )
  observed <- data.table::data.table(
    epitope = "YLQPRTFLL", n_tcr = 2622L, n_donors = 81L
  )
  out <- match_predicted_epitopes(predicted, observed)

  ## Both predicted epitopes appear; the miss carries NA rather than vanishing.
  expect_setequal(out$epitope, c("YLQPRTFLL", "KKKKKKKKK"))
  expect_true(is.na(out[epitope == "KKKKKKKKK"]$match_type))
  expect_equal(out[epitope == "YLQPRTFLL"]$n_donors, 81L)
})

test_that("one predicted epitope can match several observed entries", {
  predicted <- data.table::data.table(
    epitope = "PPSIPITVYYAVLERA", organism = "HSV-2", mhc_class = "II"
  )
  observed <- data.table::data.table(
    epitope = c("SLPITVYYA", "IPITVYYAV"),
    n_tcr = c(4L, 10L), n_donors = c(0L, 3L)
  )
  out <- match_predicted_epitopes(predicted, observed, min_core = 7L)
  expect_equal(nrow(out), 2L)
  ## Longest shared core ranks first.
  expect_gte(out$core_len[1], out$core_len[2])
})

test_that("shallow epitopes are matched, not filtered away", {
  ## The bug this guards against: comparing only against the analysable subset
  ## reports epitopes that DO have sequenced TCRs as having none.
  predicted <- data.table::data.table(
    epitope = c("YLQPRTFLL", "ELLHAPATV"),
    organism = "SARS-CoV-2", mhc_class = "I"
  )
  ## ELLHAPATV has just 2 TCRs from 1 donor - far below any analysis threshold,
  ## but a sequenced TCR nonetheless.
  observed <- data.table::data.table(
    epitope  = c("YLQPRTFLL", "ELLHAPATV"),
    n_tcr    = c(2622L, 2L),
    n_donors = c(81L, 1L)
  )
  out <- match_predicted_epitopes(predicted, observed)

  expect_equal(sum(!is.na(out$match_type)), 2L)
  expect_equal(out[epitope == "ELLHAPATV"]$n_tcr, 2L)
})

test_that("point variants are separated from the epitope itself", {
  ## Deep mutational scanning puts single-substitution variants of a popular
  ## epitope in the database. TCRs against YLQPRTFEL say something about
  ## cross-recognition, NOT that YLQPRTFLL itself was recognised.
  m <- classify_match("YLQPRTFLL", "YLQPRTFEL")
  expect_equal(m$match_class, "point variant")
  expect_equal(m$n_subs, 1L)

  m2 <- classify_match("YLQPRTFLL", "YRQPRTFLL")
  expect_equal(m2$match_class, "point variant")

  ## Beyond max_subs it is no longer the same epitope.
  m3 <- classify_match("YLQPRTFLL", "AAAAATFLL", max_subs = 2L)
  expect_false(identical(m3$match_class, "point variant"))
})

test_that("offset overlaps are classed as overlapping, not the same epitope", {
  ## Two distinct 9-mers shifted by one residue, with different HLA anchors.
  m <- classify_match("SVLNDILSR", "VLNDILSRL")
  expect_equal(m$match_class, "overlapping")
})

test_that("containment counts as the same epitope", {
  expect_equal(classify_match("TFEYVSQPFLMDLEG", "TFEYVSQPFLMDLE")$match_class,
               "same epitope")
  expect_equal(classify_match("QSAPHGVVF", "FPQSAPHGVVF")$match_class,
               "same epitope")
  expect_equal(classify_match("YLQPRTFLL", "YLQPRTFLL")$match_class,
               "same epitope")
})

test_that("strongest evidence is ranked first", {
  predicted <- data.table::data.table(
    epitope = "YLQPRTFLL", organism = "SARS-CoV-2", mhc_class = "I"
  )
  observed <- data.table::data.table(
    epitope  = c("YLQPRTFEL", "YLQPRTFLL", "QPRTFLLKY"),
    n_tcr    = c(4L, 2622L, 2L),
    n_donors = c(0L, 81L, 0L)
  )
  out <- match_predicted_epitopes(predicted, observed)
  expect_equal(out$match_class[1], "same epitope")
  expect_equal(out$matched_epitope[1], "YLQPRTFLL")
})

test_that("final ordering uses evidence strength, not the class name", {
  ## Regression: "overlapping" sorts before "same epitope" alphabetically, so
  ## ordering on the raw string silently inverted the intended ranking.
  predicted <- data.table::data.table(
    epitope = c("YLQPRTFLL", "SANNCTFEY"),
    organism = "SARS-CoV-2", mhc_class = "I"
  )
  observed <- data.table::data.table(
    epitope  = c("QPRTFLLKY", "YLQPRTFLL", "YLQPRTFEL", "SANNCTFEY"),
    n_tcr    = c(2L, 2622L, 4L, 36L),
    n_donors = c(0L, 81L, 0L, 9L)
  )
  out <- match_predicted_epitopes(predicted, observed)

  ## Within every predicted epitope, same-epitope rows must precede variants,
  ## which must precede overlaps.
  rank_of <- c("same epitope" = 1L, "point variant" = 2L, "overlapping" = 3L)
  for (ep in unique(out$epitope)) {
    r <- rank_of[out[epitope == ep]$match_class]
    expect_false(is.unsorted(r))
  }
})

test_that("unmatched epitopes do not inflate per-class counts", {
  ## Regression: `match_class == "x"` returns NA for unmatched rows, and
  ## uniqueN() counts that NA as a distinct epitope. The symptom was a summary
  ## reporting same_epitope = 1 alongside any_overlap = 0.
  res <- data.table::data.table(
    epitope     = c("AAA", "BBB", "CCC"),
    match_class = c("same epitope", NA_character_, NA_character_),
    match_type  = c("exact", NA_character_, NA_character_)
  )

  wrong <- res[, data.table::uniqueN(epitope[match_class == "same epitope"])]
  right <- res[, data.table::uniqueN(epitope[match_class %in% "same epitope"])]

  expect_equal(right, 1L)
  expect_gt(wrong, right)   # demonstrates the inflation the fix removes
})
