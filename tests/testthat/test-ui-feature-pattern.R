# -----------------------------------------------------------------------------
# The "select features by pattern" box of the Boxplots tab, on the demo
# dataset (helper-app.R): its matches fill the feature bucket, a new pattern
# replaces them, and the line under the box reports what matched.
# -----------------------------------------------------------------------------

# Type a pattern; the box is debounced, so wait for longer than that.
set_pattern <- function(app, pattern, regex = FALSE) {
  app$set_inputs(fb_pattern_regex_boxplot = regex, fb_pattern_boxplot = pattern)
  app$wait_for_idle(duration = 1500)
}

bucket <- function(app) {
  unlist(app$get_value(input = "fb_keep_boxplot")) %||% character(0)
}

pattern_info <- function(app) {
  app$get_value(output = "fb_pattern_info_boxplot")
}

test_that("a pattern fills the feature bucket with its matches", {
  app <- app_start("feature-pattern")
  app_load_demo(app)
  app_open_tab(app, "boxplotTab")

  set_pattern(app, "hspa1*")
  expect_identical(pattern_info(app), "7 names match.")
  expect_setequal(bucket(app), c("HSPA12A", "HSPA14", "HSPA1A", "HSPA1B", "HSPA1L", "HSPA13", "HSPA12B"))
  expect_identical(app_output_errors(app), character(0))

  # A new pattern replaces the old matches
  set_pattern(app, "CXCL?")
  expect_setequal(bucket(app), c("CXCL6", "CXCL5", "CXCL1", "CXCL3", "CXCL2", "CXCL8", "CXCL9"))

  set_pattern(app, "^CXCL1[0-9]$", regex = TRUE)
  expect_identical(pattern_info(app), "7 names match.")
  expect_setequal(bucket(app), c("CXCL16", "CXCL17", "CXCL10", "CXCL14", "CXCL13", "CXCL11", "CXCL12"))

  set_pattern(app, "CXCL1(", regex = TRUE)
  expect_identical(pattern_info(app), "Not a valid pattern.")
  expect_identical(bucket(app), character(0))

  set_pattern(app, "hspa1a")
  expect_identical(bucket(app), "HSPA1A")
  set_pattern(app, "")
  expect_identical(pattern_info(app), "")
  expect_identical(bucket(app), character(0))
  expect_identical(app_output_errors(app), character(0))
})
