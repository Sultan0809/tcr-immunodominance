## setup.R -------------------------------------------------------------------
## Check the R environment and install anything missing.
##
## If you are using the conda environment (environment.yml), almost everything
## is already present and this script will simply confirm that. It only falls
## back to install.packages() for what conda does not ship.
##
## Run once:  Rscript setup.R
## -----------------------------------------------------------------------------

## --- is this a conda environment? --------------------------------------------
conda_prefix <- Sys.getenv("CONDA_PREFIX")
in_conda <- nzchar(conda_prefix) &&
  startsWith(normalizePath(R.home(), mustWork = FALSE),
             normalizePath(conda_prefix, mustWork = FALSE))

if (in_conda) {
  message("[setup] conda environment detected: ", basename(conda_prefix))
  message("[setup] Prefer `conda install -c conda-forge r-<pkg>` over ",
          "install.packages() where possible.")
  message("[setup] Mixing the two can cause segfaults from toolchain mismatch.\n")
} else if (nzchar(conda_prefix)) {
  message("[setup] CONDA_PREFIX is set but R is not from that environment.")
  message("[setup] R in use: ", R.home())
  message("[setup] You may be running the system R by mistake.\n")
}

required <- c(
  "data.table", "dplyr", "tidyr", "stringr", "stringdist",
  "igraph", "jsonlite", "ggplot2", "patchwork", "scales",
  "here", "testthat"
)

## On an existing R-centric conda environment, typically only these are absent:
##     conda install -c conda-forge r-stringdist r-testthat r-ggseqlogo
## All three are leaf packages, so the solver will not touch r-base or any
## existing Bioconductor/Seurat stack.

## NOTE: immunarch is deliberately NOT a dependency. Its diversity and overlap
## metrics are implemented natively in R/diversity.R instead - they are short,
## and avoiding a CRAN compile against a conda R removes the single most
## fragile step in setting this project up.
optional <- c(
  ggseqlogo = "CDR3 motif logos (step 3)",
  pheatmap  = "V/J usage heatmaps (step 3)"
)

check <- function(pkgs) {
  vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
}

## --- required -----------------------------------------------------------------
have <- check(required)
if (all(have)) {
  message("[setup] required packages: all ", length(required), " present")
} else {
  missing <- required[!have]
  message("[setup] required packages missing: ", paste(missing, collapse = ", "))
  if (in_conda) {
    message("[setup] install them with:\n    conda install -c conda-forge ",
            paste0("r-", tolower(missing), collapse = " "))
    message("[setup] (falling back to install.packages() now; Ctrl-C to stop ",
            "and use conda instead)")
  }
  install.packages(missing, repos = "https://cloud.r-project.org")
  still <- required[!check(required)]
  if (length(still)) {
    stop("could not install: ", paste(still, collapse = ", "), call. = FALSE)
  }
}

## --- optional -----------------------------------------------------------------
message("\n[setup] optional packages:")
for (p in names(optional)) {
  status <- if (requireNamespace(p, quietly = TRUE)) "present" else "MISSING"
  message(sprintf("  %-12s %-8s  %s", p, status, optional[[p]]))
}

## --- OLGA (optional, via its command-line interface) ---------------------------
## OLGA is called as a CLI rather than through reticulate. It is an older
## package that does not support recent Python versions, so it lives in its own
## small conda environment and is invoked with system2(). This keeps the main
## environment's Python untouched and makes the call reproducible.
olga_bin <- Sys.which("olga-compute_pgen")
message("\n[setup] OLGA (olga-compute_pgen on PATH): ",
        if (nzchar(olga_bin)) olga_bin else "not found")
if (!nzchar(olga_bin)) {
  message("    Optional. Generation probability enriches step 2 but is not required;")
  message("    convergent recombination and publicness are used as proxies otherwise.")
  message("    To install it in an isolated environment:")
  message("        conda create -n olga -c conda-forge python=3.11 numpy")
  message("        conda activate olga && pip install olga && conda activate -")
}

message("\n[setup] R ", getRversion(), " at ", R.home())
message("[setup] Done. Next: Rscript tests/testthat.R")
