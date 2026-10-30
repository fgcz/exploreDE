# -----------------------------------------------------------------------------
# Utilities for testing Shiny-coupled helpers that reference reactives via
# lexical scope (e.g. bucket_add looks up `featureBucket` / `seqAnnoReactive`
# / `inputDataReactive` from its enclosing environment).
#
# with_mock_env() returns a copy of the function whose environment has been
# replaced with one containing the supplied bindings, parented onto the
# original environment so base / package lookups still work.
# -----------------------------------------------------------------------------

with_mock_env <- function(fn, ...) {
  mocks <- rlang::list2(...) # rlang::list2 supports !!! splicing
  mock_env <- new.env(parent = environment(fn))
  list2env(mocks, envir = mock_env)
  environment(fn) <- mock_env
  fn
}
