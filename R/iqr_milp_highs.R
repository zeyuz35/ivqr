### Meta -------------------------
###
### Title: Compute inverse quantile regression estimator using HiGHS
###
### Description: Solve a mixed integer linear program to compute the inverse
### quantile regression estimator using the HiGHS solver.
###
### Author: Omkar A. Katta (adapted for HiGHS)
###

### iqr_milp_highs -------------------------
#' Compute inverse quantile regression estimator using HiGHS
#'
#' Solve a mixed integer linear program to compute the inverse
#' quantile regression estimator using the HiGHS solver.
#'
#' @inheritParams iqr_milp_gurobi
#'
#' @return A named list of
#'  \enumerate{
#'    \item iqr: HiGHS model components
#'    \item params: HiGHS parameters used
#'    \item result: solution to MILP returned by HiGHS
#'    \item status: status of HiGHS's solution
#'    \item beta_X: coefficients on exogenous variables
#'    \item beta_Phi_plus: positive part of coefficients on instruments
#'    \item beta_Phi_minus: negative part of coefficients on instruments
#'    \item beta_D: coefficients on endogenous variables
#'    \item u: positive part of residuals
#'    \item v: negative part of residuals
#'    \item a: dual variable
#'    \item k: binary variable associated with u
#'    \item l: binary variable associated with v
#'    \item beta_Phi: coefficients on instruments
#'      (beta_Phi_plus - beta_Phi_minus)
#'    \item resid: residuals (u - v)
#'    \item objval: value of objective function (absolute value of beta_Phi)
#'  }
#'
#' @export
iqr_milp_highs <- function(
  Y,
  X,
  D,
  Z,
  Phi = linear_projection(D, X, Z),
  tau,
  weights = NULL,
  O_neg = NULL,
  O_pos = NULL,
  M = NULL,
  TimeLimit = 300,
  # Note: VarHintVal, VarHintPri, BranchPriority not supported by highs R package wrapper
  sparse = TRUE,
  params = list(primal_feasibility_tolerance = 1e-6, log_to_console = FALSE),
  start_method = NULL,
  start = NULL,
  fix = NULL,
  quietly = TRUE,
  show_progress = TRUE,
  ...
) {
  out <- list() # Initialize list of results to return

  send_note_if(paste("TimeLimit:", TimeLimit, "secs"), show_progress, message)

  # Get dimensions of data
  n <- length(Y)
  n_D <- nrow(D)
  n_X <- nrow(X)
  n_Z <- nrow(Z)
  n_Phi <- nrow(Phi)
  p_D <- ncol(D)
  p_X <- ncol(X)
  p_Z <- ncol(Z)
  p_Phi <- ncol(Phi)

  # Ensure that number of observations is the same
  stopifnot(all.equal(n, n_D))
  stopifnot(all.equal(n, n_X))
  stopifnot(all.equal(n, n_Z))
  stopifnot(all.equal(n, n_Phi))

  if (!is.null(weights)) {
    stopifnot(length(weights) == n)
    stopifnot(all(weights >= 0))
  }

  # If there are no endogeneous variables, return quantile regression results:
  if (p_D == 0) {
    msg <- paste("p_D is 0 -- running QR instead of IQR MILP...")
    send_note_if(msg, TRUE, warning)
    qr <- quantreg::rq(Y ~ X - 1, tau = tau, weights = weights)
    out <- qr
    return(out)
  }

  # Create vector of 1s
  ones <- rep(1, n)

  out$Phi <- Phi # by default, Phi = projection of D on X and Z

  if (is.null(M)) {
    # by default, M = 2 * max(resid from QR of Y on X and D)
    if (p_X == 0) {
      max_qr <- max(abs(
        quantreg::rq(Y ~ D - 1, tau = tau, weights = weights)$residuals
      ))
    } else if (p_D == 0) {
      max_qr <- max(abs(
        quantreg::rq(Y ~ X - 1, tau = tau, weights = weights)$residuals
      ))
    } else {
      max_qr <- max(abs(
        quantreg::rq(Y ~ X + D - 1, tau = tau, weights = weights)$residuals
      ))
    }
    M <- 2 * max_qr
  }
  out$M <- M

  num_decision_vars <- p_X + 2 * p_Phi + p_D + 5 * n

  # Objective: minimize absolute value of \beta_Phi
  obj <- c(
    rep(0, p_X), # beta_X
    rep(1, p_Phi), # beta_Phi_plus
    rep(1, p_Phi), # beta_Phi_minus
    rep(0, p_D), # beta_D
    rep(0, n), # u
    rep(0, n), # v
    rep(0, n), # a
    rep(0, n), # k
    rep(0, n)
  ) # l
  stopifnot(length(obj) == num_decision_vars)

  # Fix decision variables according to `fix`
  if (is.null(fix)) {
    A_fix <- c()
    b_fix <- c()
    sense_fix <- c()
  } else {
    stopifnot(length(fix) - num_decision_vars == 0)
    not_na <- !is.na(fix)
    A_fix <- diag(1, num_decision_vars)[not_na, ]
    b_fix <- fix[not_na]
    sense_fix <- rep('=', sum(not_na))
  }

  # Primal Feasibility Constraint (11)
  A_pf <- cbind(
    X, # beta_X
    Phi, # beta_Phi_plus
    -Phi, # beta_Phi_minus
    D, # beta_D
    diag(1, nrow = n), # u
    -diag(1, nrow = n), # v
    diag(0, nrow = n), # a
    diag(0, nrow = n), # k
    diag(0, nrow = n)
  ) # l
  b_pf <- Y
  sense_pf <- rep("=", n)

  # Dual Feasibility Constraint (13) and (14)
  if (!is.null(weights)) {
    W_X <- sweep(X, 1, weights, "*")
  } else {
    W_X <- X
  }
  A_df_X <- cbind(
    matrix(0, nrow = p_X, ncol = p_X), # beta_X
    matrix(0, nrow = p_X, ncol = p_Phi), # beta_Phi_plus
    matrix(0, nrow = p_X, ncol = p_Phi), # beta_Phi_minus
    matrix(0, nrow = p_X, ncol = p_D), # beta_D
    matrix(0, nrow = p_X, ncol = n), # u
    matrix(0, nrow = p_X, ncol = n), # v
    t(W_X), # a
    matrix(0, nrow = p_X, ncol = n), # k
    matrix(0, nrow = p_X, ncol = n)
  ) # l
  b_df_X <- (1 - tau) * t(W_X) %*% ones
  sense_df_X <- rep("=", p_X)

  if (!is.null(weights)) {
    W_Phi <- sweep(Phi, 1, weights, "*")
  } else {
    W_Phi <- Phi
  }
  A_df_Phi <- cbind(
    matrix(0, nrow = p_Phi, ncol = p_X), # beta_X
    matrix(0, nrow = p_Phi, ncol = p_Phi), # beta_Phi_plus
    matrix(0, nrow = p_Phi, ncol = p_Phi), # beta_Phi_minus
    matrix(0, nrow = p_Phi, ncol = p_D), # beta_D
    matrix(0, nrow = p_Phi, ncol = n), # u
    matrix(0, nrow = p_Phi, ncol = n), # v
    t(W_Phi), # a
    matrix(0, nrow = p_Phi, ncol = n), # k
    matrix(0, nrow = p_Phi, ncol = n)
  ) # l
  b_df_Phi <- (1 - tau) * t(W_Phi) %*% ones
  sense_df_Phi <- rep("=", p_Phi)

  # Complementary Slackness (16) and (17)
  A_cs_uk <- cbind(
    matrix(0, nrow = n, ncol = p_X), # beta_X
    matrix(0, nrow = n, ncol = p_Phi), # beta_Phi_plus
    matrix(0, nrow = n, ncol = p_Phi), # beta_Phi_minus
    matrix(0, nrow = n, ncol = p_D), # beta_D
    diag(1, nrow = n, ncol = n), # u
    matrix(0, nrow = n, ncol = n), # v
    matrix(0, nrow = n, ncol = n), # a
    -M * diag(1, nrow = n, ncol = n), # k
    matrix(0, nrow = n, ncol = n)
  ) # l
  b_cs_uk <- rep(0, n)
  sense_cs_uk <- rep("<=", n)

  A_cs_vl <- cbind(
    matrix(0, nrow = n, ncol = p_X), # beta_X
    matrix(0, nrow = n, ncol = p_Phi), # beta_Phi_plus
    matrix(0, nrow = n, ncol = p_Phi), # beta_Phi_minus
    matrix(0, nrow = n, ncol = p_D), # beta_D
    matrix(0, nrow = n, ncol = n), # u
    diag(1, nrow = n, ncol = n), # v
    matrix(0, nrow = n, ncol = n), # a
    matrix(0, nrow = n, ncol = n), # k
    -M * diag(1, nrow = n, ncol = n)
  ) # l
  b_cs_vl <- rep(0, n)
  sense_cs_vl <- rep("<=", n)

  A_cs_ak <- cbind(
    matrix(0, nrow = n, ncol = p_X), # beta_X
    matrix(0, nrow = n, ncol = p_Phi), # beta_Phi_plus
    matrix(0, nrow = n, ncol = p_Phi), # beta_Phi_minus
    matrix(0, nrow = n, ncol = p_D), # beta_D
    matrix(0, nrow = n, ncol = n), # u
    matrix(0, nrow = n, ncol = n), # v
    diag(1, nrow = n, ncol = n), # a
    -diag(1, nrow = n, ncol = n), # k
    matrix(0, nrow = n, ncol = n)
  ) # l
  b_cs_ak <- rep(0, n)
  sense_cs_ak <- rep(">=", n)

  A_cs_al <- cbind(
    matrix(0, nrow = n, ncol = p_X), # beta_X
    matrix(0, nrow = n, ncol = p_Phi), # beta_Phi_plus
    matrix(0, nrow = n, ncol = p_Phi), # beta_Phi_minus
    matrix(0, nrow = n, ncol = p_D), # beta_D
    matrix(0, nrow = n, ncol = n), # u
    matrix(0, nrow = n, ncol = n), # v
    diag(1, nrow = n, ncol = n), # a
    matrix(0, nrow = n, ncol = n), # k
    diag(1, nrow = n, ncol = n)
  ) # l
  b_cs_al <- rep(1, n)
  sense_cs_al <- rep("<=", n)

  # Pre-processing: fix residuals of outliers
  O_neg <- sort(O_neg)
  O_pos <- sort(O_pos)
  O <- c(O_neg, O_pos) # indices of fixed residuals
  if (!is.null(O)) {
    fixed <- rep(0, n)
    fixed[O] <- 1
    fixed_mat <- diag(fixed)

    A_pp_a <- cbind(
      matrix(0, nrow = n, ncol = p_X), # beta_X
      matrix(0, nrow = n, ncol = p_Phi), # beta_Phi_plus
      matrix(0, nrow = n, ncol = p_Phi), # beta_Phi_minus
      matrix(0, nrow = n, ncol = p_D), # beta_D
      matrix(0, nrow = n, ncol = n), # u
      matrix(0, nrow = n, ncol = n), # v
      fixed_mat, # a
      matrix(0, nrow = n, ncol = n), # k
      matrix(0, nrow = n, ncol = n)
    ) # l
    b_pp_a <- rep(0, n)
    b_pp_a[O_pos] <- 1
    b_pp_a[O_neg] <- 0
    sense_pp_a <- rep("=", n)

    A_pp_k <- cbind(
      matrix(0, nrow = n, ncol = p_X), # beta_X
      matrix(0, nrow = n, ncol = p_Phi), # beta_Phi_plus
      matrix(0, nrow = n, ncol = p_Phi), # beta_Phi_minus
      matrix(0, nrow = n, ncol = p_D), # beta_D
      matrix(0, nrow = n, ncol = n), # u
      matrix(0, nrow = n, ncol = n), # v
      matrix(0, nrow = n, ncol = n), # a
      fixed_mat, # k
      matrix(0, nrow = n, ncol = n)
    ) # l
    b_pp_k <- rep(0, n)
    b_pp_k[O_pos] <- 1
    b_pp_k[O_neg] <- 0
    sense_pp_k <- rep("=", n)

    A_pp_l <- cbind(
      matrix(0, nrow = n, ncol = p_X), # beta_X
      matrix(0, nrow = n, ncol = p_Phi), # beta_Phi_plus
      matrix(0, nrow = n, ncol = p_Phi), # beta_Phi_minus
      matrix(0, nrow = n, ncol = p_D), # beta_D
      matrix(0, nrow = n, ncol = n), # u
      matrix(0, nrow = n, ncol = n), # v
      matrix(0, nrow = n, ncol = n), # a
      matrix(0, nrow = n, ncol = n), # k
      fixed_mat
    ) # l
    b_pp_l <- rep(0, n)
    b_pp_l[O_pos] <- 0
    b_pp_l[O_neg] <- 1
    sense_pp_l <- rep("=", n)
  } else {
    A_pp_a <- A_pp_k <- A_pp_l <- b_pp_a <- b_pp_k <- b_pp_l <- sense_pp_a <- sense_pp_k <- sense_pp_l <- c()
  }

  # Combine constraints
  A_all <- rbind(
    A_fix,
    A_pf,
    A_df_X,
    A_df_Phi,
    A_cs_uk,
    A_cs_vl,
    A_cs_ak,
    A_cs_al,
    A_pp_a,
    A_pp_k,
    A_pp_l
  )
  rhs_all <- c(
    b_fix,
    b_pf,
    b_df_X,
    b_df_Phi,
    b_cs_uk,
    b_cs_vl,
    b_cs_ak,
    b_cs_al,
    b_pp_a,
    b_pp_k,
    b_pp_l
  )
  sense_all <- c(
    sense_fix,
    sense_pf,
    sense_df_X,
    sense_df_Phi,
    sense_cs_uk,
    sense_cs_vl,
    sense_cs_ak,
    sense_cs_al,
    sense_pp_a,
    sense_pp_k,
    sense_pp_l
  )

  if (sparse) {
    A_all <- methods::as(A_all, "sparseMatrix")
  }

  # Map sense to HiGHS lhs/rhs
  lhs_highs <- rep(-Inf, length(sense_all))
  rhs_highs <- rep(Inf, length(sense_all))

  lhs_highs[sense_all == "="] <- rhs_all[sense_all == "="]
  rhs_highs[sense_all == "="] <- rhs_all[sense_all == "="]

  rhs_highs[sense_all == "<="] <- rhs_all[sense_all == "<="]

  lhs_highs[sense_all == ">="] <- rhs_all[sense_all == ">="]

  # Bounds
  lb <- c(
    rep(-Inf, p_X),
    rep(0, 2 * p_Phi),
    rep(-Inf, p_D),
    rep(0, 2 * n),
    rep(0, 3 * n)
  )
  ub <- c(
    rep(Inf, p_X),
    rep(Inf, 2 * p_Phi),
    rep(Inf, p_D),
    rep(Inf, 2 * n),
    rep(1, 3 * n)
  )

  # Types
  vtype <- c(rep("C", p_X + 2 * p_Phi + p_D + 3 * n), rep("I", 2 * n))
  types_highs <- ifelse(vtype == "C", 1L, 2L)

  # Solve
  params$time_limit <- TimeLimit
  result <- highs::highs_solve(
    L = obj,
    lower = lb,
    upper = ub,
    A = A_all,
    lhs = lhs_highs,
    rhs = rhs_highs,
    types = types_highs,
    control = params
  )

  status_map <- list(
    "Optimal" = "OPTIMAL",
    "Infeasible" = "INFEASIBLE",
    "Unbounded" = "UNBOUNDED",
    "Time limit" = "TIME_LIMIT"
  )
  status <- status_map[[result$status_message]]
  if (is.null(status)) {
    status <- result$status_message
  }

  out$iqr <- list(
    obj = obj,
    A = A_all,
    lhs = lhs_highs,
    rhs = rhs_highs,
    lb = lb,
    ub = ub,
    types = types_highs
  )
  out$params <- params
  out$result <- result
  out$status <- status

  if (status %in% c("OPTIMAL", "SUBOPTIMAL", "Optimal")) {
    answer <- result$primal_solution
    if (p_X > 0) {
      out$beta_X <- answer[1:p_X]
    } else {
      out$beta_X <- NA
    }
    if (p_Phi > 0) {
      out$beta_Phi_plus <- answer[(p_X + 1):(p_X + p_Phi)]
      out$beta_Phi_minus <- answer[(p_X + p_Phi + 1):(p_X + 2 * p_Phi)]
    } else {
      out$beta_Phi_plus <- NA
      out$beta_Phi_minus <- NA
    }
    if (p_D > 0) {
      out$beta_D <- answer[(p_X + 2 * p_Phi + 1):(p_X + 2 * p_Phi + p_D)]
    } else {
      out$beta_D <- NA
    }
    out$u <- answer[(p_X + 2 * p_Phi + p_D + 1):(p_X + 2 * p_Phi + p_D + n)]
    out$v <- answer[
      (p_X + 2 * p_Phi + p_D + n + 1):(p_X + 2 * p_Phi + p_D + 2 * n)
    ]
    out$a <- answer[
      (p_X + 2 * p_Phi + p_D + 2 * n + 1):(p_X + 2 * p_Phi + p_D + 3 * n)
    ]
    out$k <- answer[
      (p_X + 2 * p_Phi + p_D + 3 * n + 1):(p_X + 2 * p_Phi + p_D + 4 * n)
    ]
    out$l <- answer[
      (p_X + 2 * p_Phi + p_D + 4 * n + 1):(p_X + 2 * p_Phi + p_D + 5 * n)
    ]

    out$beta_Phi <- out$beta_Phi_plus - out$beta_Phi_minus
    out$resid <- out$u - out$v
    out$objval <- result$objective_value
  }

  return(out)
}
