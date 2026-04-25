library(testthat)
library(ivqr)

test_that("iqr_milp works with highs solver", {
  if (!requireNamespace("highs", quietly = TRUE)) {
    skip("HiGHS not installed")
  }

  set.seed(1)
  n <- 50
  X <- matrix(1, n, 1)
  Z <- matrix(rnorm(n), n, 1)
  D <- 0.5 * Z + rnorm(n)
  Y <- 1 + D + rnorm(n)

  # Solve with HiGHS
  fit_highs <- iqr_milp(
    Y,
    X,
    D,
    Z,
    tau = 0.5,
    solver = "highs",
    quietly = TRUE,
    show_progress = FALSE
  )

  expect_type(fit_highs$beta_D, "double")
  expect_true(!is.na(fit_highs$beta_D))

  # If Gurobi is also available, compare results
  if (requireNamespace("gurobi", quietly = TRUE)) {
    fit_gurobi <- iqr_milp(
      Y,
      X,
      D,
      Z,
      tau = 0.5,
      solver = "gurobi",
      quietly = TRUE,
      show_progress = FALSE
    )
    expect_equal(fit_highs$beta_D, fit_gurobi$beta_D, tolerance = 1e-4)
  }
})

test_that("preprocess_iqr_milp works with highs solver", {
  if (!requireNamespace("highs", quietly = TRUE)) {
    skip("HiGHS not installed")
  }

  set.seed(1)
  n <- 100
  X <- matrix(1, n, 1)
  Z <- matrix(rnorm(n), n, 1)
  D <- 0.5 * Z + rnorm(n)
  Y <- 1 + D + rnorm(n)

  fit_pre_highs <- preprocess_iqr_milp(
    Y,
    D,
    X,
    Z,
    tau = 0.5,
    solver = "highs",
    show_iterations = FALSE
  )

  expect_type(fit_pre_highs$final_fit$beta_D, "double")
})

test_that("test_stat works with highs solver", {
  if (!requireNamespace("highs", quietly = TRUE)) {
    skip("HiGHS not installed")
  }

  set.seed(1)
  n <- 100
  X <- matrix(1, n, 1)
  Z <- matrix(rnorm(n), n, 1)
  D <- 0.5 * Z + rnorm(n)
  Y <- 1 + D + rnorm(n)

  ts_highs <- test_stat(
    beta_D_null = 1,
    beta_X_null = NA,
    Y = Y,
    X = X,
    D = D,
    Z = Z,
    tau = 0.5,
    solver = "highs",
    show_progress = FALSE,
    print_results = FALSE
  )

  expect_type(ts_highs$test_stat, "double")
})
