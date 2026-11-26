# -----------------------------------------------------------------------------
# ?data= links (helper-app.R): a link the app cannot read shows a message and
# keeps the session.
# -----------------------------------------------------------------------------

test_that("a proteomics ?data= link to a file that cannot be read shows a message", {
  app <- app_session("broken-proteomics-link")
  # The link points at pStore; the file does not exist (or pStore is not
  # reachable), so reading it fails either way.
  app$run_js("window.location.search = '?data=Proteomics/exploreDE-test/does-not-exist.rds';")
  app$wait_for_js("document.querySelectorAll('.modal-dialog').length > 0", timeout = 60e3)
  expect_match(app$get_js("document.querySelector('.modal-dialog').innerText"), "can't be opened")
  expect_true(app$get_js("Shiny.shinyapp.isConnected()"))
})

test_that("an .h5ad ?data= link that is not on gStore shows a message", {
  app <- app_session("missing-h5ad-link")
  # Not on gStore (or the course gStore), so the .h5ad route falls through to
  # the "does not exist" message rather than trying to read it.
  app$run_js("window.location.search = '?data=p99999/exploreDE-test/does-not-exist.h5ad';")
  app$wait_for_js("document.querySelectorAll('.modal-dialog').length > 0", timeout = 60e3)
  expect_match(app$get_js("document.querySelector('.modal-dialog').innerText"), "does not exist")
  expect_true(app$get_js("Shiny.shinyapp.isConnected()"))
})
