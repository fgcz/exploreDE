# KEGG and GO helpers: organism codes, pathway lookups, GO annotation columns.

# Helper functions

#' kegg_species_label
#' Maps a `ref_build`/`refBuild` string's leading "Genus_species" segment
#' (e.g. `"Homo_sapiens/GENCODE/..."`) to the KEGG-organism display name the
#' pathway tabs use to decide whether KEGG lookups are supported. Self-
#' contained so exploreDE does not need `ezRun::getSpecies()` just for this —
#' `ezRun` is otherwise optional (`Suggests`), and a fresh install without it
#' would hard-error on the four call sites this replaces (app.R, server-kegg.R,
#' server-goTile.R, server-heatmap2.R). Only the three organisms those call
#' sites ever recognised are mapped; anything else returns `NULL`, matching
#' their existing `switch(spp, "Human"=, "Mouse"=, "Rat"=, NULL)` behaviour.
#'
#' @param ref_build A `ref_build`/`refBuild` string, or `NA`/`NULL`.
#' @return `"Human"`, `"Mouse"`, `"Rat"`, or `NULL`.
#' @noRd
kegg_species_label <- function(ref_build) {
  if (is.null(ref_build) || is.na(ref_build) || !nzchar(ref_build)) {
    return(NULL)
  }
  genus_species <- strsplit(ref_build, "/", fixed = TRUE)[[1]][1]
  switch(genus_species, Homo_sapiens = "Human", Mus_musculus = "Mouse", Rattus_norvegicus = "Rat", NULL)
}

#' kegg_installed
#' The KEGG views need pathview and KEGGREST, both in Suggests. Every KEGG
#' entry point (menu items, the KEGG inputs of the tile and heatmap tabs) is
#' shown only when this is `TRUE`. Installed, not loaded: the sidebar menu
#' asks on every render, and loading pathview takes about a second, so it
#' loads when a KEGG view first uses it.
#' @return `TRUE` when both packages are installed.
#' @noRd
kegg_installed <- function() {
  nzchar(system.file(package = "pathview")) && nzchar(system.file(package = "KEGGREST"))
}

#' kegg_org_code
#' KEGG organism code of a reference build, from `kegg_species_label()`.
#' @param ref_build A `ref_build`/`refBuild` string, or `NA`/`NULL`.
#' @return `"hsa"`, `"mmu"`, `"rno"`, or `NULL` for any other species.
#' @noRd
kegg_org_code <- function(ref_build) {
  switch(kegg_species_label(ref_build) %||% "", Human = "hsa", Mouse = "mmu", Rat = "rno", NULL)
}

#' kegg_fetch_pathway
#' Name and member genes of a KEGG pathway, from the KEGG REST server. The
#' only network access of the KEGG tile and heatmap views; its errors (no
#' connection, unknown pathway) are for the caller to report.
#' @param pathway_id Pathway id with organism prefix, e.g. `"hsa04110"`.
#' @param org KEGG organism code, e.g. `"hsa"`.
#' @return `list(name, entrez)`: the pathway name without its organism suffix,
#'   and the Entrez ids of its genes.
#' @noRd
kegg_fetch_pathway <- function(pathway_id, org) {
  info <- KEGGREST::keggGet(pathway_id)
  links <- KEGGREST::keggLink(org, paste0("path:", pathway_id))
  list(name = sub(" - .*$", "", info[[1]]$NAME), entrez = gsub(paste0(org, ":"), "", links))
}

#' kegg_map_genes
#' Feature ids of a DE table whose genes are in a KEGG pathway. Matches the
#' table's `entrez_id` column when it has one; otherwise maps the Entrez ids
#' to gene symbols and Ensembl ids with `pathview::eg2id()` and matches
#' `gene_name` and `gene_id`.
#' @param de A DE table with `feature_id`.
#' @param entrez Entrez ids, as from `kegg_fetch_pathway()`.
#' @param org KEGG organism code.
#' @return Character vector of feature ids.
#' @noRd
kegg_map_genes <- function(de, entrez, org) {
  if ("entrez_id" %in% colnames(de)) {
    return(de$feature_id[de$entrez_id %in% entrez])
  }
  ids <- pathview::eg2id(entrez, category = c("SYMBOL", "ENSEMBL"), org = org)
  unique(c(
    de$feature_id[de$gene_name %in% stats::na.omit(ids[, "SYMBOL"])],
    de$feature_id[de$gene_id %in% stats::na.omit(ids[, "ENSEMBL"])]
  ))
}

#' go_annotation_columns
#' The GO annotation columns a feature table carries, out of `GO BP`, `GO MF`,
#' and `GO CC`. Any subset is valid: Non-reference organisms often have only
#' one ontology. An empty result means the GO tabs are hidden, as KEGG is for
#' unsupported species.
#'
#' @param df A feature or DE table (`data.frame`).
#' @return Character vector of the GO columns present, possibly empty.
#' @noRd
go_annotation_columns <- function(df) {
  intersect(c("GO BP", "GO MF", "GO CC"), colnames(df))
}

#' features_in_go_term
#' Which rows of a feature table are annotated to one GO term, looking in
#' whichever GO columns the table has (see `go_annotation_columns()`). Each
#' column holds `"; "`-separated GO ids; `NA` means no annotation.
#'
#' @param df A feature or DE table (`data.frame`).
#' @param go_id A GO id, e.g. `"GO:0006955"`.
#' @return Logical vector, one per row; all `FALSE` without GO columns.
#' @noRd
features_in_go_term <- function(df, go_id) {
  hit <- rep(FALSE, nrow(df))
  for (col in go_annotation_columns(df)) {
    x <- as.character(df[[col]])
    hit <- hit | (!is.na(x) & grepl(go_id, x, fixed = TRUE))
  }
  hit
}
