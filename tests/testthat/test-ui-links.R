# -----------------------------------------------------------------------------
# ?data= links (helper-app.R): a link the app cannot read shows a message and
# keeps the session.
# -----------------------------------------------------------------------------

test_that("a proteomics ?data= link to a file that cannot be read shows a message", {
  app <- app_start("broken-proteomics-link")
  # The link points at pStore; the file does not exist (or pStore is not
  # reachable), so reading it fails either way.
  app$run_js("window.location.search = '?data=Proteomics/exploreDE-test/does-not-exist.rds';")
  app$wait_for_js("document.querySelectorAll('.modal-dialog').length > 0", timeout = 60e3)
  expect_match(app$get_js("document.querySelector('.modal-dialog').innerText"), "can't be opened")
  expect_true(app$get_js("Shiny.shinyapp.isConnected()"))
})
