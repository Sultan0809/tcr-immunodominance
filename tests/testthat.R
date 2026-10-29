## Run the test suite:  Rscript tests/testthat.R
library(testthat)
library(data.table)

source(file.path("R", "utils.R"))
source(file.path("R", "load_vdjdb.R"))
source(file.path("R", "load_vdjdb_chunks.R"))
source(file.path("R", "annotate_antigen.R"))
source(file.path("R", "diversity.R"))
source(file.path("R", "features.R"))
source(file.path("R", "cdr3_network.R"))
source(file.path("R", "publicness_null.R"))
source(file.path("R", "match_epitopes.R"))

test_dir("tests/testthat")
