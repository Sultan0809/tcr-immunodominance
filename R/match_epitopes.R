## R/match_epitopes.R --------------------------------------------------------
## Match externally predicted epitopes against VDJdb's observed repertoires.
##
## WHY NOT EXACT MATCHING
##
## Predicted epitopes and database entries rarely align exactly. A class II
## 15-mer and a database 14-mer describing the same binding core are the same
## epitope for TCR purposes; so are a 9-mer and the 16-mer that contains it.
## Matching is therefore done on the longest common substring, with its length
## reported so the user can judge each call rather than trusting a flag.
##
## This answers a narrow question: has a human TCR ever been sequenced against
## this epitope, or against something that overlaps it? A miss means no data,
## not non-immunogenicity.
## -----------------------------------------------------------------------------

#' Longest common substring length between two strings
#'
#' Classic dynamic programme. Sequences here are tens of residues, so the
#' quadratic cost is irrelevant.
#'
#' @param a,b Character scalars.
#' @return A list with `len` and the matched substring `seq`.
lcs_length <- function(a, b) {
  if (is.na(a) || is.na(b) || !nzchar(a) || !nzchar(b)) {
    return(list(len = 0L, seq = ""))
  }
  x <- strsplit(a, "")[[1]]; y <- strsplit(b, "")[[1]]
  n <- length(x); m <- length(y)
  dp <- matrix(0L, n + 1L, m + 1L)
  best <- 0L; end_i <- 0L

  for (i in seq_len(n)) {
    for (j in seq_len(m)) {
      if (x[i] == y[j]) {
        dp[i + 1L, j + 1L] <- dp[i, j] + 1L
        if (dp[i + 1L, j + 1L] > best) {
          best <- dp[i + 1L, j + 1L]; end_i <- i
        }
      }
    }
  }
  list(len = best,
       seq = if (best > 0L) paste(x[(end_i - best + 1L):end_i], collapse = "") else "")
}

#' Classify how a predicted epitope relates to a database epitope
#'
#' @param pred,db Character scalars.
#' @param min_core Minimum shared-core length to count as a match.
#' @return A list with `type`, `core`, `core_len`.
classify_match <- function(pred, db, min_core = 7L, max_subs = 2L) {
  mk <- function(type, cls, core, len, subs = NA_integer_) {
    list(type = type, match_class = cls, core = core, core_len = len,
         n_subs = subs)
  }

  ## --- same epitope ---------------------------------------------------------
  ## The database entry and the prediction describe one epitope: identical, or
  ## one contains the other (a 9-mer inside a 15-mer, a 14-mer inside a 15-mer).
  if (identical(pred, db)) {
    return(mk("exact", "same epitope", pred, nchar(pred), 0L))
  }
  if (grepl(db, pred, fixed = TRUE)) {
    return(mk("db within predicted", "same epitope", db, nchar(db)))
  }
  if (grepl(pred, db, fixed = TRUE)) {
    return(mk("predicted within db", "same epitope", pred, nchar(pred)))
  }

  ## --- point variant --------------------------------------------------------
  ## Same length, a few substitutions: a mutational-scanning entry or a strain
  ## variant of the SAME epitope. TCRs here speak to cross-recognition and
  ## escape tolerance, not to recognition of the predicted sequence itself.
  if (nchar(pred) == nchar(db)) {
    subs <- sum(strsplit(pred, "")[[1]] != strsplit(db, "")[[1]])
    if (subs <= max_subs) {
      l <- lcs_length(pred, db)
      return(mk(sprintf("point variant (%d sub)", subs), "point variant",
                l$seq, l$len, as.integer(subs)))
    }
  }

  ## --- overlapping ----------------------------------------------------------
  ## A shared stretch, but offset: often a DIFFERENT epitope from the same
  ## region, with its own HLA restriction. Weakest evidence - report it, but
  ## never count it as recognition of the predicted epitope.
  l <- lcs_length(pred, db)
  if (l$len >= min_core) {
    return(mk("shared core", "overlapping", l$seq, l$len))
  }

  mk(NA_character_, NA_character_, NA_character_, l$len)
}

#' Match a set of predicted epitopes against observed VDJdb repertoires
#'
#' @param predicted A `data.table` with at least an `epitope` column.
#' @param observed A per-epitope feature table (from `epitope_features()` or a
#'   downstream script), keyed on `epitope`.
#' @param min_core Minimum shared-core length.
#' @return One row per (predicted, matched) pair, plus unmatched predicted
#'   epitopes with `NA` match fields, so nothing silently disappears.
match_predicted_epitopes <- function(predicted, observed, min_core = 7L) {
  stopifnot("epitope" %in% names(predicted), "epitope" %in% names(observed))

  obs_eps <- unique(observed$epitope)
  log_msg("matching %d predicted epitopes against %d observed (min core = %d)",
          nrow(predicted), length(obs_eps), min_core)

  rows <- lapply(seq_len(nrow(predicted)), function(i) {
    p <- predicted[i]
    hits <- lapply(obs_eps, function(o) {
      m <- classify_match(p$epitope, o, min_core = min_core)
      if (is.na(m$type)) return(NULL)
      data.table::data.table(
        matched_epitope = o,
        match_type      = m$type,
        match_class     = m$match_class,
        shared_core     = m$core,
        core_len        = m$core_len,
        n_subs          = m$n_subs
      )
    })
    hits <- data.table::rbindlist(Filter(Negate(is.null), hits))

    if (!nrow(hits)) {
      return(cbind(p, data.table::data.table(
        matched_epitope = NA_character_, match_type = NA_character_,
        match_class = NA_character_, shared_core = NA_character_,
        core_len = NA_integer_, n_subs = NA_integer_
      )))
    }
    ## Strongest evidence first: same epitope, then point variant, then
    ## overlapping; within a class, the longest shared stretch.
    hits[, class_rank := data.table::fcase(
      match_class == "same epitope",  1L,
      match_class == "point variant", 2L,
      default = 3L
    )]
    data.table::setorder(hits, class_rank, -core_len)
    hits[, class_rank := NULL]
    cbind(p[rep(1L, nrow(hits))], hits)
  })

  out <- data.table::rbindlist(rows, fill = TRUE)

  ## Attach the repertoire features of whatever was matched.
  obs_cols <- intersect(
    c("epitope", "protein", "mhc_a", "n_tcr", "n_cdr3", "n_donors", "n_studies",
      "public_fraction", "cluster_public_fraction", "public_rare_mean",
      "public_z", "clonality", "v_entropy"),
    names(observed)
  )
  obs_sub <- unique(observed[, ..obs_cols])
  data.table::setnames(obs_sub, "epitope", "matched_epitope")

  out <- merge(out, obs_sub, by = "matched_epitope", all.x = TRUE, sort = FALSE)

  n_same <- data.table::uniqueN(out[match_class == "same epitope"]$epitope)
  n_var  <- data.table::uniqueN(out[match_class == "point variant"]$epitope)
  n_any  <- data.table::uniqueN(out[!is.na(match_type)]$epitope)
  n_pred <- data.table::uniqueN(predicted$epitope)
  log_msg("  same epitope: %d/%d | point variant: %d | any overlap: %d",
          n_same, n_pred, n_var, n_any)

  data.table::setcolorder(out, c("epitope", "organism", "mhc_class",
                                 "matched_epitope", "match_class", "match_type",
                                 "shared_core", "core_len", "n_subs"))

  ## Order by evidence strength, NOT by the class name. Sorting the string
  ## alphabetically puts "overlapping" ahead of "same epitope", which is the
  ## opposite of the intended ranking.
  out[, .class_rank := data.table::fcase(
    match_class == "same epitope",  1L,
    match_class == "point variant", 2L,
    match_class == "overlapping",   3L,
    default = 4L
  )]
  data.table::setorder(out, organism, .class_rank, -core_len, epitope)
  out[, .class_rank := NULL]
  out[]
}
