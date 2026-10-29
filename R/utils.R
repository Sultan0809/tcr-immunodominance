## R/utils.R -----------------------------------------------------------------
## Small shared helpers: project paths, logging, safe IO.
## Sourced by scripts; has no side effects beyond defining functions.
## -----------------------------------------------------------------------------

#' Project-root-relative path
#'
#' Thin wrapper so the pipeline works from RStudio, `Rscript`, or a Quarto
#' render without hard-coded absolute paths.
#'
#' @param ... Path components, as in `file.path()`.
#' @return An absolute path (character scalar).
proj_path <- function(...) {
  if (requireNamespace("here", quietly = TRUE)) {
    here::here(...)
  } else {
    file.path(normalizePath(".", mustWork = FALSE), ...)
  }
}

#' Directory shortcuts
dir_raw       <- function(...) proj_path("data", "raw", ...)
dir_processed <- function(...) proj_path("data", "processed", ...)
dir_figures   <- function(...) proj_path("figures", ...)

#' Ensure a directory exists, returning its path invisibly
ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

#' Timestamped console logging
#'
#' @param ... Passed to `sprintf()` if more than one argument, else pasted.
#' @param level One of "info", "warn", "ok".
log_msg <- function(..., level = c("info", "warn", "ok")) {
  level <- match.arg(level)
  args <- list(...)
  txt <- if (length(args) > 1L) do.call(sprintf, args) else as.character(args[[1]])
  tag <- switch(level, info = "[info]", warn = "[warn]", ok = "[ ok ]")
  message(sprintf("%s %s %s", format(Sys.time(), "%H:%M:%S"), tag, txt))
  invisible(txt)
}

#' Download a file with retries, returning TRUE on success
#'
#' Does not stop on failure — the caller decides whether a failed source is
#' fatal or just one of several fallbacks.
#'
#' @param url Source URL.
#' @param dest Destination path.
#' @param tries Number of attempts.
#' @param quiet Suppress download progress.
#' @return Logical scalar: did a non-empty file land at `dest`?
try_download <- function(url, dest, tries = 2L, quiet = TRUE) {
  ensure_dir(dirname(dest))
  for (i in seq_len(tries)) {
    ok <- tryCatch({
      utils::download.file(url, dest, mode = "wb", quiet = quiet)
      TRUE
    }, error = function(e) {
      log_msg("download attempt %d/%d failed: %s", i, tries,
              conditionMessage(e), level = "warn")
      FALSE
    }, warning = function(w) {
      log_msg("download attempt %d/%d warning: %s", i, tries,
              conditionMessage(w), level = "warn")
      FALSE
    })
    if (ok && file.exists(dest) && file.info(dest)$size > 0) return(TRUE)
  }
  if (file.exists(dest) && file.info(dest)$size == 0) unlink(dest)
  FALSE
}

#' Normalise VDJdb-style column names to snake_case
#'
#' VDJdb uses dotted names (`antigen.epitope`, `vdjdb.score`) and web exports
#' sometimes use spaces or different capitalisation. This makes downstream code
#' independent of which export the user supplied.
#'
#' @param x Character vector of column names.
#' @return Character vector, lowercase snake_case.
normalise_names <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

#' Write a processed table and log it
write_processed <- function(x, filename) {
  ensure_dir(dir_processed())
  path <- dir_processed(filename)
  data.table::fwrite(x, path, sep = "\t")
  log_msg("wrote %s (%d rows x %d cols)", filename, nrow(x), ncol(x), level = "ok")
  invisible(path)
}
