## Tests for the VDJdb chunk-format loader. Offline: chunks are synthesised.

make_chunk_dir <- function() {
  d <- file.path(tempdir(), paste0("chunks_", as.integer(runif(1, 1e6, 9e6))))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

write_chunk <- function(dir, name, dt) {
  data.table::fwrite(dt, file.path(dir, name), sep = "\t")
}

test_that("load_vdjdb_chunks reshapes paired chains into long form", {
  d <- make_chunk_dir()
  write_chunk(d, "PMID_1.tsv", data.table::data.table(
    cdr3_alpha      = c("CAVRDSNYQLIW", ""),
    v_alpha         = c("TRAV3", ""),
    j_alpha         = c("TRAJ33", ""),
    cdr3_beta       = c("CASSLGQAYEQYF", "CASSIRSSYEQYF"),
    v_beta          = c("TRBV7-9", "TRBV19"),
    j_beta          = c("TRBJ2-7", "TRBJ2-7"),
    d_beta          = c("TRBD1", ""),
    species         = c("HomoSapiens", "HomoSapiens"),
    mhc_a           = c("HLA-A*02:01", "HLA-A*02:01"),
    mhc_class       = c("MHCI", "MHCI"),
    antigen_epitope = c("GILGFVFTL", "GILGFVFTL"),
    antigen_gene    = c("M", "M"),
    antigen_species = c("InfluenzaA", "InfluenzaA"),
    reference_id    = c("PMID:1", "PMID:1"),
    meta_study_id   = c("PMID:1", "PMID:1"),
    meta_subject_id = c("D1", "D2"),
    method_identification = c("tetramer-sort", "culture")
  ))

  beta <- load_vdjdb_chunks(d, chain = "TRB")
  expect_equal(nrow(beta), 2L)
  expect_true(all(beta$gene == "TRB"))
  expect_setequal(beta$cdr3, c("CASSLGQAYEQYF", "CASSIRSSYEQYF"))

  ## Only one row has an alpha chain; the empty one must be dropped.
  alpha <- load_vdjdb_chunks(d, chain = "TRA")
  expect_equal(nrow(alpha), 1L)
  expect_equal(alpha$cdr3, "CAVRDSNYQLIW")

  unlink(d, recursive = TRUE)
})

test_that("donor_key combines study and subject", {
  d <- make_chunk_dir()
  ## Same subject label in two different studies must not collapse into one donor.
  write_chunk(d, "PMID_1.tsv", data.table::data.table(
    cdr3_beta = "CASSLGQAYEQYF", v_beta = "TRBV7-9", j_beta = "TRBJ2-7",
    species = "HomoSapiens", antigen_epitope = "GILGFVFTL",
    reference_id = "PMID:1", meta_study_id = "PMID:1", meta_subject_id = "1"
  ))
  write_chunk(d, "PMID_2.tsv", data.table::data.table(
    cdr3_beta = "CASSIRSSYEQYF", v_beta = "TRBV19", j_beta = "TRBJ2-7",
    species = "HomoSapiens", antigen_epitope = "GILGFVFTL",
    reference_id = "PMID:2", meta_study_id = "PMID:2", meta_subject_id = "1"
  ))

  res <- load_vdjdb_chunks(d, chain = "TRB")
  expect_equal(nrow(res), 2L)
  expect_equal(data.table::uniqueN(res$donor_key), 2L)
  expect_true(all(res$has_subject_id))

  unlink(d, recursive = TRUE)
})

test_that("missing subject IDs yield NA donor_key, not a fake donor", {
  d <- make_chunk_dir()
  write_chunk(d, "PMID_1.tsv", data.table::data.table(
    cdr3_beta = c("CASSLGQAYEQYF", "CASSIRSSYEQYF"),
    v_beta = c("TRBV7-9", "TRBV19"), j_beta = c("TRBJ2-7", "TRBJ2-7"),
    species = "HomoSapiens", antigen_epitope = "GILGFVFTL",
    reference_id = "PMID:1", meta_subject_id = c("", NA)
  ))
  res <- load_vdjdb_chunks(d, chain = "TRB")
  ## Assert on the parsed value itself, so a failure reports what actually
  ## came back rather than just TRUE/FALSE.
  expect_equal(res$subject_id, c(NA_character_, NA_character_))
  expect_equal(res$donor_key, c(NA_character_, NA_character_))
  expect_false(any(res$has_subject_id))
  unlink(d, recursive = TRUE)
})

test_that("blank_to_na covers the spellings submitters actually use", {
  expect_equal(blank_to_na(c("D1", "", " ", NA, "NA", "N/A", "null", ".", "0")),
               c("D1", NA, NA, NA, NA, NA, NA, NA, "0"))
})

test_that("blank_to_na strips quoting left behind by quote-disabled reading", {
  ## Chunks are read with quote = "", so a quoted empty field arrives as the
  ## literal two-character string: \"\" . It must still count as missing.
  expect_equal(blank_to_na('""'), NA_character_)
  expect_equal(blank_to_na(c('"D1"', "'D2'", '"  "', '""')),
               c("D1", "D2", NA, NA))
  ## An internal quote is not a wrapper and must survive.
  expect_equal(blank_to_na('say "hi"'), 'say "hi"')
})

test_that("chunks with differing columns bind without error", {
  d <- make_chunk_dir()
  write_chunk(d, "PMID_1.tsv", data.table::data.table(
    cdr3_beta = "CASSLGQAYEQYF", v_beta = "TRBV7-9",
    species = "HomoSapiens", antigen_epitope = "GILGFVFTL",
    reference_id = "PMID:1"
  ))
  ## Second chunk has extra and missing fields relative to the first.
  write_chunk(d, "PMID_2.tsv", data.table::data.table(
    cdr3_beta = "CASSIRSSYEQYF", v_beta = "TRBV19", j_beta = "TRBJ2-7",
    species = "HomoSapiens", antigen_epitope = "NLVPMVATV",
    reference_id = "PMID:2", meta_tissue = "PBMC", comment = "x"
  ))
  res <- load_vdjdb_chunks(d, chain = "TRB")
  expect_equal(nrow(res), 2L)
  expect_true(all(c("tissue", "j_segm") %in% names(res)))
  unlink(d, recursive = TRUE)
})

test_that("non-amino-acid CDR3s are dropped", {
  d <- make_chunk_dir()
  write_chunk(d, "PMID_1.tsv", data.table::data.table(
    cdr3_beta = c("CASSLGQAYEQYF", "CASS_BAD*SEQ", "casslower"),
    v_beta = "TRBV7-9", species = "HomoSapiens",
    antigen_epitope = "GILGFVFTL", reference_id = "PMID:1"
  ))
  res <- load_vdjdb_chunks(d, chain = "TRB")
  expect_equal(nrow(res), 1L)
  expect_equal(res$cdr3, "CASSLGQAYEQYF")
  unlink(d, recursive = TRUE)
})

test_that("chunk_confidence flags direct binding evidence only", {
  dt <- data.table::data.table(
    method_identification = c("tetramer-sort", "multimer-sort", "culture",
                              "antigen-loaded-targets", NA, "structure")
  )
  got <- chunk_confidence(dt)
  expect_equal(got, c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE))
})

test_that("find_chunks_dir prefers the directory with most TSVs", {
  root <- file.path(tempdir(), paste0("root_", as.integer(runif(1, 1e6, 9e6))))
  big <- file.path(root, "a", "chunks")
  small <- file.path(root, "b", "chunks")
  dir.create(big, recursive = TRUE); dir.create(small, recursive = TRUE)
  for (i in 1:3) file.create(file.path(big, sprintf("c%d.tsv", i)))
  file.create(file.path(small, "c1.tsv"))

  expect_equal(normalizePath(find_chunks_dir(root)), normalizePath(big))
  unlink(root, recursive = TRUE)
})

test_that("blank_to_na treats empty and whitespace as missing", {
  expect_equal(blank_to_na(c("D1", "", "  ", NA, "0")),
               c("D1", NA, NA, NA, "0"))
  expect_equal(blank_to_na(character(0)), character(0))
})

test_that("is_natural_repertoire excludes in vitro display", {
  dt <- data.table::data.table(
    method_identification = c("tetramer-sort",
                              "phage display, magnetic beads",
                              "yeast display",
                              "dextramer-sort",
                              NA)
  )
  expect_equal(is_natural_repertoire(dt), c(TRUE, FALSE, FALSE, TRUE, TRUE))
})

test_that("display-derived records do not inflate donor counts", {
  d <- make_chunk_dir()
  write_chunk(d, "PMID_display.tsv", data.table::data.table(
    cdr3_beta = c("CASSLGQAYEQYF", "CASSIRSSYEQYF"),
    v_beta = "TRBV7-9", species = "HomoSapiens",
    antigen_epitope = "GILGFVFTL", reference_id = "PMID:9",
    method_identification = c("phage display, magnetic beads", "tetramer-sort")
  ))
  res <- load_vdjdb_chunks(d, chain = "TRB")
  expect_equal(sum(is_natural_repertoire(res)), 1L)
  unlink(d, recursive = TRUE)
})
