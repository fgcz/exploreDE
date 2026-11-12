test_that("detect_read_count_column prefers exact names, then 'Read' + 'Count', then 'Read'", {
  expect_identical(detect_read_count_column(c("Name", "Read Count", "Reads")), "Reads")
  expect_identical(detect_read_count_column(c("Name", "PF Read Count", "Read File")), "PF Read Count")
  expect_identical(detect_read_count_column(c("Name", "Read Mapped", "Raw Read Pairs")), "Raw Read Pairs")
  expect_null(detect_read_count_column(c("Name", "Condition")))
})

test_that("library_size_text gives the range of library sizes in millions", {
  ds <- data.frame(Name = c("a", "b", "c"), Reads = c(12.34e6, 20e6, 31.06e6))
  expect_identical(
    library_size_text(ds),
    "Individual library size ranged from 12.3 million to 31.1 million reads. "
  )
  expect_identical(library_size_text(data.frame(Name = "a", Condition = "x")), "")
  expect_identical(library_size_text(data.frame(Name = "a", Reads = "n/a")), "")
})

test_that("build_methods_docx writes every section and adds GO only when present", {
  skip_if_not_installed("officer")
  sections <- list(data_processing = "dp", trimming = "tr", mapping = "mp", de_analysis = "de", go_analysis = NULL)
  text <- function(doc) officer::docx_summary(doc)$text
  doc <- build_methods_docx(sections, list(a = "Citation A", b = "Citation B"), date = as.Date("2026-11-12"))
  expect_true(all(
    c("Methods", "Data Processing", "Differential Expression Analysis", "References", "Citation B") %in% text(doc)
  ))
  expect_false("Gene Ontology Enrichment Analysis" %in% text(doc))
  expect_true(any(grepl("12 November 2026", text(doc))))
  sections$go_analysis <- "go"
  expect_true("Gene Ontology Enrichment Analysis" %in% text(build_methods_docx(sections, list())))
})
