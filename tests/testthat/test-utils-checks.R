test_that("lm_check_flag() accepts only a single TRUE or FALSE", {
  expect_silent(lm_check_flag(TRUE, "progress"))
  expect_silent(lm_check_flag(FALSE, "progress"))

  expect_error(lm_check_flag(NA, "progress"), "TRUE.*FALSE")
  expect_error(lm_check_flag(c(TRUE, FALSE), "progress"), "TRUE.*FALSE")
  expect_error(lm_check_flag(1, "progress"), "TRUE.*FALSE")
  expect_error(lm_check_flag("yes", "progress"), "TRUE.*FALSE")
})

test_that("lm_progress_arg() returns the label or FALSE for purrr::map()", {
  expect_identical(lm_progress_arg(TRUE, "Fitting mixed models"), "Fitting mixed models")
  expect_false(lm_progress_arg(FALSE, "Fitting mixed models"))
})

test_that("the progress option falls back to interactive() and can be overridden", {
  old <- options(tidylandmark.progress = NULL)

  expect_identical(
    getOption("tidylandmark.progress", interactive()),
    interactive()
  )

  options(tidylandmark.progress = FALSE)
  expect_false(getOption("tidylandmark.progress", interactive()))

  options(old)
})
