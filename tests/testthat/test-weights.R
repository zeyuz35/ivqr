library(testthat)
library(ivqr)

test_that("iqr_milp works with identity weights", {
  # Skip if gurobi is not installed
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    skip("Gurobi not installed")
  }

  set.seed(1)
  n <- 50
  X <- matrix(1, n, 1)
  Z <- matrix(rnorm(n), n, 1)
  D <- 0.5 * Z + rnorm(n)
  Y <- 1 + D + rnorm(n)

  # Standard call
  fit1 <- iqr_milp(Y, X, D, Z, tau = 0.5, quietly = TRUE, show_progress = FALSE)

  # Identity weights
  fit2 <- iqr_milp(
    Y,
    X,
    D,
    Z,
    tau = 0.5,
    weights = rep(1, n),
    quietly = TRUE,
    show_progress = FALSE
  )

  expect_equal(fit1$beta_D, fit2$beta_D, tolerance = 1e-5)
  expect_equal(fit1$beta_X, fit2$beta_X, tolerance = 1e-5)
})

test_that("iqr_milp works with random weights", {
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    skip("Gurobi not installed")
  }

  set.seed(1)
  n <- 50
  X <- matrix(1, n, 1)
  Z <- matrix(rnorm(n), n, 1)
  D <- 0.5 * Z + rnorm(n)
  Y <- 1 + D + rnorm(n)
  w <- runif(n, 0.5, 1.5)

  fit_w <- iqr_milp(
    Y,
    X,
    D,
    Z,
    tau = 0.5,
    weights = w,
    quietly = TRUE,
    show_progress = FALSE
  )

  expect_type(fit_w$beta_D, "double")
  expect_true(!is.na(fit_w$beta_D))
})

test_that("test_stat works with weights", {
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    skip("Gurobi not installed")
  }

  set.seed(1)
  n <- 100
  X <- matrix(1, n, 1)
  Z <- matrix(rnorm(n), n, 1)
  D <- 0.5 * Z + rnorm(n)
  Y <- 1 + D + rnorm(n)
  w <- runif(n, 0.5, 1.5)

  # Test stat under null beta_D = 1
  ts_w <- test_stat(
    beta_D_null = 1,
    beta_X_null = NA,
    Y = Y,
    X = X,
    D = D,
    Z = Z,
    tau = 0.5,
    weights = w,
    show_progress = FALSE,
    print_results = FALSE
  )

  expect_type(ts_w$test_stat, "double")
  expect_type(ts_w$p_val, "double")
})
