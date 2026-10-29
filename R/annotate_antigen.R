## R/annotate_antigen.R ------------------------------------------------------
## Map epitope peptides onto reference proteins by sequence, not by label.
##
## VDJdb's `antigen.gene` field is free text and inconsistent ("Spike", "S",
## "surface glycoprotein"). Locating each epitope *inside* the actual reference
## protein sequence is both more reliable and gives the residue position, which
## step 3 needs. `antigen.gene` is kept only as a fallback.
## -----------------------------------------------------------------------------

## Reference proteins fetched from UniProt.
## `expect` is a regex checked against the FASTA header so a wrong or
## repurposed accession is caught loudly instead of silently mislabelling data.
.REFERENCE_PROTEINS <- data.table::data.table(
  organism  = c(rep("SARS-CoV-2", 6L), "HSV-2"),
  protein   = c("S", "N", "M", "E", "ORF3a", "ORF1ab", "gD"),
  structural = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE),
  accession = c("P0DTC2", "P0DTC9", "P0DTC5", "P0DTC4",
                "P0DTC3", "P0DTD1", "P03172"),
  expect    = c("Spike", "Nucleoprotein", "Membrane", "Envelope small",
                "ORF3a|Protein 3a", "Replicase", "glycoprotein D|Envelope glycoprotein D")
)

#' Fetch one protein sequence from UniProt
#'
#' @param accession UniProt accession, e.g. "P0DTC2".
#' @param cache_dir Where to cache the FASTA.
#' @return A list with `header` and `sequence`, or `NULL` on failure.
fetch_uniprot <- function(accession, cache_dir = dir_raw("uniprot")) {
  ensure_dir(cache_dir)
  dest <- file.path(cache_dir, paste0(accession, ".fasta"))

  if (!file.exists(dest) || file.info(dest)$size == 0) {
    url <- sprintf("https://rest.uniprot.org/uniprotkb/%s.fasta", accession)
    if (!try_download(url, dest)) {
      log_msg("could not fetch UniProt %s", accession, level = "warn")
      return(NULL)
    }
  }

  lines <- readLines(dest, warn = FALSE)
  if (!length(lines) || !startsWith(lines[1], ">")) {
    log_msg("UniProt %s: malformed FASTA", accession, level = "warn")
    return(NULL)
  }
  list(
    header   = sub("^>", "", lines[1]),
    sequence = toupper(gsub("[^A-Za-z]", "", paste(lines[-1], collapse = "")))
  )
}

#' Build the reference proteome table
#'
#' @param proteins Subset of `.REFERENCE_PROTEINS` to fetch. Defaults to all.
#' @param verify Check each FASTA header against the expected protein name.
#' @return A `data.table` with one row per successfully fetched protein.
build_reference_proteome <- function(proteins = .REFERENCE_PROTEINS, verify = TRUE) {
  rows <- lapply(seq_len(nrow(proteins)), function(i) {
    p <- proteins[i]
    fa <- fetch_uniprot(p$accession)
    if (is.null(fa)) return(NULL)

    if (verify && !grepl(p$expect, fa$header, ignore.case = TRUE)) {
      log_msg(paste0("accession %s header does not match expected '%s'.\n",
                     "         header: %s\n",
                     "         Check the accession in .REFERENCE_PROTEINS before trusting ",
                     "any annotation for %s."),
              p$accession, p$expect, fa$header, p$protein, level = "warn")
    }

    data.table::data.table(
      organism   = p$organism,
      protein    = p$protein,
      structural = p$structural,
      accession  = p$accession,
      length     = nchar(fa$sequence),
      header     = fa$header,
      sequence   = fa$sequence
    )
  })

  out <- data.table::rbindlist(Filter(Negate(is.null), rows))
  if (!nrow(out)) stop("No reference proteins could be fetched.", call. = FALSE)
  log_msg("reference proteome: %d proteins (%s)",
          nrow(out), paste(out$protein, collapse = ", "), level = "ok")
  out[]
}

#' Locate epitopes within reference proteins
#'
#' Exact substring match. Epitopes from variant strains will not match the
#' reference and are returned with `protein = NA`; the count is logged so the
#' miss rate is visible rather than silent.
#'
#' @param epitopes Character vector of peptide sequences (duplicates fine).
#' @param proteome Output of `build_reference_proteome()`.
#' @return A `data.table`: one row per unique epitope, with `protein`,
#'   `organism`, `structural`, `start`, `end`, `rel_position`.
annotate_epitope_position <- function(epitopes, proteome) {
  uniq <- unique(epitopes[!is.na(epitopes) & nzchar(epitopes)])
  if (!length(uniq)) {
    return(data.table::data.table(
      epitope = character(0), organism = character(0), protein = character(0),
      structural = logical(0), start = integer(0), end = integer(0),
      rel_position = numeric(0)
    ))
  }

  hits <- lapply(uniq, function(ep) {
    pos <- vapply(proteome$sequence, function(s) {
      as.integer(regexpr(ep, s, fixed = TRUE))
    }, integer(1))
    idx <- which(pos > 0L)

    if (!length(idx)) {
      return(data.table::data.table(
        epitope = ep, organism = NA_character_, protein = NA_character_,
        structural = NA, start = NA_integer_, end = NA_integer_,
        rel_position = NA_real_
      ))
    }
    ## An epitope matching more than one protein is almost always a short
    ## peptide inside a polyprotein; keep the shortest protein as the most
    ## specific assignment.
    i <- idx[which.min(proteome$length[idx])]
    data.table::data.table(
      epitope      = ep,
      organism     = proteome$organism[i],
      protein      = proteome$protein[i],
      structural   = proteome$structural[i],
      start        = pos[i],
      end          = pos[i] + nchar(ep) - 1L,
      rel_position = (pos[i] + nchar(ep) / 2) / proteome$length[i]
    )
  })

  out <- data.table::rbindlist(hits)
  n_hit <- sum(!is.na(out$protein))
  log_msg("epitope mapping: %d/%d matched a reference protein (%.1f%%)",
          n_hit, nrow(out), 100 * n_hit / nrow(out))
  if (n_hit < nrow(out)) {
    log_msg("unmatched epitopes are from variant strains or other antigens; ",
            level = "info")
    log_msg("  they keep VDJdb's `epitope_gene` label as a fallback.", level = "info")
  }
  out[]
}

#' Attach epitope annotation to a VDJdb table
#'
#' @param vdjdb Output of `load_vdjdb()`.
#' @param proteome Output of `build_reference_proteome()`.
#' @return `vdjdb` with annotation columns joined on `epitope`.
annotate_vdjdb <- function(vdjdb, proteome) {
  ann <- annotate_epitope_position(vdjdb$epitope, proteome)
  out <- merge(vdjdb, ann, by = "epitope", all.x = TRUE, sort = FALSE)

  ## Fallback label when no sequence match was possible.
  out[, protein_label := data.table::fcoalesce(protein, epitope_gene, "unassigned")]
  out[]
}
