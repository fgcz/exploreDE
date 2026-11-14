# -----------------------------------------------------------------------------
# Colours / Labels tab (helper-app.R): the colour template download lists one
# row per group level with its colour and label, and importing it back
# applies the colours.
# -----------------------------------------------------------------------------

test_that("the colour template round-trips through download and import", {
  app <- app_session("colour-template")
  app_load_demo(app)
  app_open_tab(app, "colourTab")
  app$wait_for_js("!!$('#dlColourTemplate').attr('href')")

  template <- openxlsx::read.xlsx(app$get_download("dlColourTemplate"))
  expect_identical(colnames(template), c("Level", "Colour", "Label"))
  expect_gt(nrow(template), 0)
  expect_true(all(grepl("^GroupColour", template$Level)))

  template$Colour[1] <- "#123456"
  f <- file.path(withr::local_tempdir(), "colours.xlsx")
  openxlsx::write.xlsx(template, f)
  app$upload_file(importColourFile = f)
  app$wait_for_idle()
  expect_identical(tolower(app$get_value(input = template$Level[1])), "#123456")
  expect_identical(app_output_errors(app), character(0))
})
