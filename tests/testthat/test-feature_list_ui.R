test_that("feature_list_ui renders feature ids as text, not HTML", {
  html <- as.character(feature_list_ui(c("<img src=x onerror=alert(1)>", "TP53")))
  expect_false(grepl("<img", html, fixed = TRUE))
  expect_true(grepl("&lt;img src=x onerror=alert(1)&gt;", html, fixed = TRUE))
  expect_true(grepl("<div>TP53</div>", html, fixed = TRUE))
})
