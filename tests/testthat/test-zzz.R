context("test-zzz")

test_that("get and set global variables strings", {
  expect_true(is.na(get_global_variable('regression_test_variable')))
  expect_message(get_global_variable('regression_test_variable'))

  expect_silent(set_global_variable('regression_test_variable', 5))
  expect_equal(get_global_variable('regression_test_variable'), 5)
} )

