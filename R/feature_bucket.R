# Feature bucket and feature selection: bucket updates, pasted lists, name patterns, symbol lookup.

#' bucket_add
#' The feature bucket after adding features: `current`, then those of `new`
#' that are `known` and not yet in it, in their order.
#' @param current Feature ids in the bucket.
#' @param new Feature ids to add (for CRISPR guides collapsed to genes, gene
#'   names).
#' @param known Ids that may enter the bucket.
#' @return Character vector.
#' @noRd
bucket_add <- function(current, new, known) {
  c(current, setdiff(new[new %in% known & new != ""], current))
}

#' bucket_remove
#' The feature bucket without `ids`.
#' @param current Feature ids in the bucket.
#' @param ids Feature ids to remove.
#' @return Character vector.
#' @noRd
bucket_remove <- function(current, ids) {
  current[!current %in% ids]
}

#' parse_gene_text
#' Split a textarea string of gene names by whitespace / comma / slash / newline / tab.
#' @param text Character scalar
#' @return Character vector of non-empty tokens
#' @noRd
parse_gene_text <- function(text) {
  if (is.null(text) || text == "") {
    return(character(0))
  }
  genes <- strsplit(text, "\\ +|\\/+|\\,+|\n+|\r+|\n\r+|\t+")[[1]]
  genes[genes != ""]
}

#' match_feature_pattern
#' The candidates a "select features by pattern" box matches, ignoring case.
#' By default the pattern is a name with the wildcards `*` (any characters)
#' and `?` (one character): `Hsp*` is every name starting with Hsp, a bare
#' name only that name, and every other character is literal, so protein IDs
#' like `sp|P02768|ALBU_HUMAN` match as typed. With `regex = TRUE` it is an
#' extended regular expression (`^Hsp[0-9]`, `rik$`, `Gm1|Gm2`; no Perl
#' lookarounds), unanchored. An empty
#' or invalid pattern matches nothing; an invalid one is flagged by
#' `n_matched = NA`.
#' @param pattern Character scalar, the box's value
#' @param candidates Character vector of names to match against
#' @param regex Take the pattern as a regular expression
#' @return The matches, in `candidates` order, with attribute `n_matched`
#' @noRd
match_feature_pattern <- function(pattern, candidates, regex = FALSE) {
  pattern <- trimws(pattern %||% "")
  if (!nzchar(pattern) || !length(candidates)) {
    return(structure(character(0), n_matched = 0L))
  }
  if (!isTRUE(regex)) {
    pattern <- gsub("([.|()\\[\\]{}+^$\\\\])", "\\\\\\1", pattern, perl = TRUE)
    pattern <- gsub("?", ".", gsub("*", ".*", pattern, fixed = TRUE), fixed = TRUE)
    pattern <- paste0("^", pattern, "$")
  }
  # R's default (TRE) engine, not PCRE: it does not backtrack, so a pattern
  # such as (.+)+X$ cannot hold the R process (and every session on it) for
  # minutes. An invalid pattern makes TRE warn, then error: either means no
  # match can be computed.
  hit <- tryCatch(
    grepl(pattern, candidates, ignore.case = TRUE),
    warning = function(w) NULL,
    error = function(e) NULL
  )
  if (is.null(hit)) {
    return(structure(character(0), n_matched = NA_integer_))
  }
  matched <- unique(candidates[hit])
  structure(matched, n_matched = length(matched))
}

#' cap_pattern_features
#' The features a pattern box adds: its matched names resolved to feature ids
#' (a protein ID or gene symbol can stand for many peptides or sites), at most
#' `max_features` of them.
#' @param matched Names from `match_feature_pattern()`
#' @param resolve Function turning names into feature ids
#' @param max_features Keep at most this many feature ids
#' @return The feature ids, with attribute `n_features`: how many before the cap
#' @noRd
cap_pattern_features <- function(matched, resolve, max_features = 1000L) {
  if (!length(matched)) {
    return(structure(character(0), n_features = 0L))
  }
  ids <- resolve(as.character(matched))
  structure(utils::head(ids, max_features), n_features = length(ids))
}

#' fb_pattern_message
#' The line shown under a pattern box: nothing for an empty box, otherwise how
#' many names matched and whether the feature cap cut them.
#' @param pattern Character scalar, the box's value
#' @param matched The value of `match_feature_pattern()` for it
#' @param features The value of `cap_pattern_features()` for those matches
#' @return Character scalar
#' @noRd
fb_pattern_message <- function(pattern, matched, features = NULL) {
  if (!nzchar(trimws(pattern %||% ""))) {
    return("")
  }
  n <- attr(matched, "n_matched") %||% length(matched)
  if (is.na(n)) {
    return("Not a valid pattern.")
  }
  n_features <- attr(features, "n_features") %||% length(features)
  if (n_features > length(features)) {
    return(sprintf(
      "%d names match (%d features); the first %d features were added.",
      n,
      n_features,
      length(features)
    ))
  }
  if (n == 1) "1 name matches." else sprintf("%d names match.", n)
}

#' text_input_on_blur
#' A text box that sends its value when it loses focus or on Enter, not on
#' every keystroke, where the installed shiny supports `updateOn`.
#' @param ... Passed to `shiny::textInput()`
#' @return A shiny tag
#' @noRd
text_input_on_blur <- function(...) {
  args <- list(...)
  if ("updateOn" %in% names(formals(shiny::textInput))) {
    args$updateOn <- "blur"
  }
  do.call(shiny::textInput, args)
}

#' fb_pattern_ui
#' The "select features by pattern" controls of a feature-bucket tab: the
#' pattern box, the "Regular expression" tick box and the line reporting the
#' matches.
#' @param tab The tab's feature-bucket id: "boxplot", "volcano" or "heatmap"
#' @return A shiny tag list
#' @noRd
fb_pattern_ui <- function(tab) {
  shiny::tagList(
    text_input_on_blur(
      inputId = paste0("fb_pattern_", tab),
      label = "Or, select features by pattern:",
      placeholder = "Hsp* or ^Hsp[0-9] (regular expression)"
    ),
    shiny::checkboxInput(
      inputId = paste0("fb_pattern_regex_", tab),
      label = "Regular expression (otherwise * and ? are wildcards)",
      value = FALSE
    ),
    shiny::textOutput(paste0("fb_pattern_info_", tab))
  )
}

#' resolve_features
#' Feature ids for the tokens a user entered in "Search by" mode `mode`: gene
#' symbols or protein ids expanded through their lookup, plus tokens that are
#' themselves known feature ids. In feature_id mode, or without the lookup,
#' the tokens as entered.
#' @param tokens Character vector of user input tokens.
#' @param mode One of "feature_id", "gene_name", "protein_Id".
#' @param lookups List with `gene` and `protein` lookups (name -> feature ids),
#'   as `eds$proteomics_choices(kind)$lookup`.
#' @param known Feature ids of the DE table.
#' @return Character vector of feature ids.
#' @noRd
resolve_features <- function(tokens, mode = "feature_id", lookups = list(), known = character()) {
  lookup <- switch(mode, gene_name = lookups$gene, protein_Id = lookups$protein, NULL)
  if (is.null(lookup)) {
    return(tokens)
  }
  unique(c(unlist(lookup[tokens], use.names = FALSE), tokens[tokens %in% known]))
}
