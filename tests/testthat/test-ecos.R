library(testthat)
library(ivqr)

test_that("iqr_milp works with ecos solver", {
  if (!requireNamespace("ECOSolveR", quietly = TRUE)) {
    skip("ECOSolveR not installed")
  }

  set.seed(1)
  n <- 10
  X <- matrix(1, n, 1)
  Z <- matrix(rnorm(n), n, 1)
  D <- 0.5 * Z + rnorm(n)
  Y <- 1 + D + rnorm(n)

  # Solve with ECOS
  fit_ecos <- iqr_milp(
    Y,
    X,
    D,
    Z,
    tau = 0.5,
    solver = "ecos",
    quietly = TRUE,
    show_progress = FALSE
  )

  expect_type(fit_ecos$beta_D, "double")
  expect_true(!is.na(fit_ecos$beta_D))

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
    expect_equal(fit_ecos$beta_D, fit_gurobi$beta_D, tolerance = 1e-4)
  }
})

test_that("miqcp_proj works with ecos solver", {
  if (!requireNamespace("ECOSolveR", quietly = TRUE)) {
    skip("ECOSolveR not installed")
  }

  set.seed(42)
  n <- 5
  data <- chen_lee(n = n, p_D = 1)
  Y <- data$Y
  D <- data$D
  X <- data$X
  Z <- data$Z
  tau <- 0.5

  # We need an initial fit to get residuals
  fit <- iqr_milp(
    Y,
    X,
    D,
    Z,
    tau = tau,
    solver = "ecos",
    quietly = TRUE,
    show_progress = FALSE
  )

  # Subvector projection
  res_ecos <- miqcp_proj(
    projection_index = 1,
    endogeneous = TRUE,
    alpha = 0.1,
    sense = "min",
    Y = Y,
    X = X,
    D = D,
    Z = Z,
    tau = tau,
    residuals = fit$resid,
    solver = "ecos",
    show_progress = FALSE,
    quietly = TRUE
  )

  expect_type(res_ecos$objval, "double")

  if (requireNamespace("gurobi", quietly = TRUE)) {
    res_gurobi <- miqcp_proj(
      projection_index = 1,
      endogeneous = TRUE,
      alpha = 0.1,
      sense = "min",
      Y = Y,
      X = X,
      D = D,
      Z = Z,
      tau = tau,
      residuals = fit$resid,
      solver = "gurobi",
      show_progress = FALSE,
      quietly = TRUE
    )
    expect_equal(res_ecos$objval, res_gurobi$objval, tolerance = 1e-4)
  }
})
