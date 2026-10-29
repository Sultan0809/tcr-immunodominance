## Tests for the pure parsing/annotation functions.
## These run offline - no VDJdb or UniProt download required.
## Run with:  testthat::test_dir("tests/testthat")

test_that("normalise_names handles VDJdb and web-export spellings", {
  expect_equal(normalise_names("antigen.epitope"), "antigen_epitope")
  expect_equal(normalise_names("MHC A"), "mhc_a")
  expect_equal(normalise_names("  Score  "), "score")
  expect_equal(normalise_names("V"), "v")
  expect_equal(normalise_names(c("CDR3", "complex.id")), c("cdr3", "complex_id"))
})

test_that("extract_meta_field preserves vector alignment", {
  meta <- c(
    '{"study.id":"PMID:1","subject.id":"D7","tissue":"PBMC"}',
    '{"study.id":"PMID:2","subject.id":"","tissue":""}',
    NA_character_,
    '{"study.id":"PMID:3","subject.id":"D9"}'
  )
  got <- extract_meta_field(meta, "subject.id")

  ## Critical: output length must equal input length, with NA in the right slots.
  expect_length(got, 4L)
  expect_equal(got, c("D7", NA, NA, "D9"))

  expect_equal(extract_meta_field(meta, "study.id"),
               c("PMID:1", "PMID:2", NA, "PMID:3"))
})

test_that("extract_meta_field returns all-NA for an absent field", {
  meta <- c('{"study.id":"PMID:1"}', '{"study.id":"PMID:2"}')
  expect_true(all(is.na(extract_meta_field(meta, "donor.MHC"))))
})

test_that("extract_meta_field copes with empty input", {
  expect_length(extract_meta_field(character(0), "subject.id"), 0L)
})

test_that("annotate_epitope_position locates peptides and reports misses", {
  fake_proteome <- data.table::data.table(
    organism   = c("TestOrg", "TestOrg"),
    protein    = c("ShortProt", "LongProt"),
    structural = c(TRUE, FALSE),
    accession  = c("X1", "X2"),
    length     = c(20L, 40L),
    header     = c("h1", "h2"),
    sequence   = c("MKLAYTGWQNPVRSDEFGHI",
                   "AAAAAMKLAYTGWQNPVRSDEFGHICCCCCDDDDDEEEEE")
  )

  res <- annotate_epitope_position(c("YTGWQ", "NOTPRESENT"), fake_proteome)

  expect_equal(nrow(res), 2L)

  hit <- res[epitope == "YTGWQ"]
  ## Present in both proteins; the shorter one wins as the more specific call.
  expect_equal(hit$protein, "ShortProt")
  expect_equal(hit$start, 5L)
  expect_equal(hit$end, 9L)
  expect_true(hit$structural)

  miss <- res[epitope == "NOTPRESENT"]
  expect_true(is.na(miss$protein))
  expect_true(is.na(miss$start))
})

test_that("annotate_epitope_position handles empty and NA input", {
  fake_proteome <- data.table::data.table(
    organism = "TestOrg", protein = "P", structural = TRUE, accession = "X",
    length = 10L, header = "h", sequence = "ACDEFGHIKL"
  )
  expect_equal(nrow(annotate_epitope_position(character(0), fake_proteome)), 0L)
  expect_equal(nrow(annotate_epitope_position(c(NA, ""), fake_proteome)), 0L)
})

test_that("is_vdjdb_table validates by header, not filename", {
  tmp <- tempfile(fileext = ".tsv")

  ## A real VDJdb header.
  writeLines("complex.id\tgene\tcdr3\tv.segm\tantigen.epitope\tvdjdb.score", tmp)
  expect_true(is_vdjdb_table(tmp))

  ## A decoy with a VDJdb-looking name but unrelated columns - this is the
  ## file that broke the first version of find_vdjdb_file().
  writeLines("name\tlogin\tavatar\tcontributions", tmp)
  expect_false(is_vdjdb_table(tmp))

  ## Web-export spelling must also validate.
  writeLines("Gene\tCDR3\tV\tJ\tEpitope\tScore", tmp)
  expect_true(is_vdjdb_table(tmp))

  ## Empty file.
  writeLines(character(0), tmp)
  expect_false(is_vdjdb_table(tmp))

  unlink(tmp)
})

test_that("is_vdjdb_table returns FALSE for a missing file", {
  expect_false(is_vdjdb_table(file.path(tempdir(), "definitely-not-here.tsv")))
})
