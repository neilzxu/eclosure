test_that("timing benchmark: closedeBH on 1e6 e-values", {
  skip_if(
    Sys.getenv("ECLOSURE_RUN_TIMING_TESTS") != "1",
    "Set ECLOSURE_RUN_TIMING_TESTS=1 to run timing benchmarks"
  )

  set.seed(1)
  e <- rexp(1e6, rate = 1 / 5)

  timing <- system.time({
    res <- closedeBH(e)
  })

  expect_true(is.numeric(res) || is.integer(res))
  expect_length(res, 1)
  expect_gte(res, 0)
  message(sprintf(
    "closedeBH timing (1e6 e-values): %.3f seconds",
    timing[["elapsed"]]
  ))
})


test_that("timing benchmark: closedBY on 1e5 p-values", {
  skip_if(
    Sys.getenv("ECLOSURE_RUN_TIMING_TESTS") != "1",
    "Set ECLOSURE_RUN_TIMING_TESTS=1 to run timing benchmarks"
  )

  set.seed(1)
  p <- runif(1e5)^4

  timing <- system.time({
    res <- closedBY(p)
  })

  expect_true(is.numeric(res) || is.integer(res))
  expect_length(res, 1)
  expect_gte(res, 0)
  message(sprintf(
    "closedBY timing (1e5 p-values): %.3f seconds",
    timing[["elapsed"]]
  ))
})
