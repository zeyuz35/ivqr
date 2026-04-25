### Meta -------------------------
###
### Title: Compute inverse quantile regression estimator using ECOS
###
### Description: Solve a mixed integer linear program to compute the inverse
### quantile regression estimator using the ECOS solver.
###
### Author: Omkar A. Katta (adapted for ECOS)
###

### iqr_milp_ecos -------------------------
#' Compute inverse quantile regression estimator using ECOS
#'
#' Solve a mixed integer linear program to compute the inverse
#' quantile regression estimator using the ECOS solver.
#'
#' @inheritParams iqr_milp_gurobi
#'
#' @return A named list of
#'  \enumerate{
#'    \item iqr: ECOS model components
#'    \item params: ECOS parameters used
#'    \item result: solution to MILP returned by ECOS
#'    \item status: status of ECOS's solution
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
iqr_milp_ecos <- function(
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
  sparse = TRUE,
  params = ECOSolveR::ecos.control(),
  start_method = NULL,
  start = NULL,
  fix = NULL,
  quietly = TRUE,
  show_progress = TRUE,
  ...
) {
  out <- list() # Initialize list of results to return

  if (missing(params) || is.null(params)) {
    params <- ECOSolveR::ecos.control(mi_max_iters = 10000L)
  }
  # Ensure verbose is set according to quietly
  params$verbose <- ifelse(quietly, 0L, 1L)

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
    A_fix <- diag(1, num_decision_vars)[not_na, , drop = FALSE]
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

  stopifnot(ncol(A_pf) == num_decision_vars)
  stopifnot(nrow(A_pf) == n)
  stopifnot(length(b_pf) == n)
  stopifnot(length(sense_pf) == n)
  msg <- paste("Primal Feasibility Complete.")
  send_note_if(msg, show_progress, message)

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

  stopifnot(ncol(A_df_X) == num_decision_vars)
  stopifnot(nrow(A_df_X) == p_X)
  stopifnot(length(b_df_X) == p_X)
  stopifnot(length(sense_df_X) == p_X)
  msg <- paste("Dual Feasibility for X Complete.")
  send_note_if(msg, show_progress, message)

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

  stopifnot(ncol(A_df_Phi) == num_decision_vars)
  stopifnot(nrow(A_df_Phi) == p_Phi)
  stopifnot(length(b_df_Phi) == p_Phi)
  stopifnot(length(sense_df_Phi) == p_Phi)
  msg <- paste("Dual Feasibility for Phi Complete.")
  send_note_if(msg, show_progress, message)

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

  stopifnot(ncol(A_cs_uk) == num_decision_vars)
  stopifnot(nrow(A_cs_uk) == n)
  stopifnot(length(b_cs_uk) == n)
  stopifnot(length(sense_cs_uk) == n)
  msg <- paste("Complementary Slackness for u and k Complete.")
  send_note_if(msg, show_progress, message)

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

  stopifnot(ncol(A_cs_vl) == num_decision_vars)
  stopifnot(nrow(A_cs_vl) == n)
  stopifnot(length(b_cs_vl) == n)
  stopifnot(length(sense_cs_vl) == n)
  msg <- paste("Complementary Slackness for v and l Complete.")
  send_note_if(msg, show_progress, message)

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

  stopifnot(ncol(A_cs_ak) == num_decision_vars)
  stopifnot(nrow(A_cs_ak) == n)
  stopifnot(length(b_cs_ak) == n)
  stopifnot(length(sense_cs_ak) == n)
  msg <- paste("Complementary Slackness for a and k Complete.")
  send_note_if(msg, show_progress, message)

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

  stopifnot(ncol(A_cs_al) == num_decision_vars)
  stopifnot(nrow(A_cs_al) == n)
  stopifnot(length(b_cs_al) == n)
  stopifnot(length(sense_cs_al) == n)
  msg <- paste("Complementary Slackness for a and l Complete.")
  send_note_if(msg, show_progress, message)

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

  # Bounds
  # ECOSolveR takes bounds as inequalities Gx <= h
  # x >= lb => -x <= -lb
  # x <= ub => x <= ub
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

  # Convert lb and ub to inequalities
  # For finite lb and ub
  I_mat <- diag(1, num_decision_vars)

  idx_lb <- which(lb != -Inf)
  G_lb <- -I_mat[idx_lb, , drop = FALSE]
  h_lb <- -lb[idx_lb]

  idx_ub <- which(ub != Inf)
  G_ub <- I_mat[idx_ub, , drop = FALSE]
  h_ub <- ub[idx_ub]

  # Map sense to ECOSolveR format
  # Equality: A_ecos %*% x = b_ecos
  A_ecos <- A_all[sense_all == "=", , drop = FALSE]
  b_ecos <- rhs_all[sense_all == "="]

  # Inequality: G_ecos %*% x + s = h_ecos, s in K
  G_ecos_orig <- rbind(
    A_all[sense_all == "<=", , drop = FALSE],
    -A_all[sense_all == ">=", , drop = FALSE]
  )
  h_ecos_orig <- c(
    rhs_all[sense_all == "<="],
    -rhs_all[sense_all == ">="]
  )

  G_ecos <- rbind(G_ecos_orig, G_lb, G_ub)
  h_ecos <- c(h_ecos_orig, h_lb, h_ub)

  if (sparse) {
    A_ecos <- methods::as(A_ecos, "sparseMatrix")
    G_ecos <- methods::as(G_ecos, "sparseMatrix")
  }

  # Boolean variables: k and l are the last 2n
  bool_vars <- (num_decision_vars - 2 * n + 1):num_decision_vars

  # Solve
  result <- ECOSolveR::ECOS_csolve(
    c = obj,
    G = G_ecos,
    h = h_ecos,
    dims = list(l = length(h_ecos), q = NULL, e = 0L),
    A = A_ecos,
    b = b_ecos,
    bool_vars = bool_vars,
    control = params
  )

  if (!quietly) {
    print(result$infostring)
    print(result$retcodes)
  }

  status_map <- list(
    "0" = "OPTIMAL",
    "1" = "INFEASIBLE",
    "2" = "UNBOUNDED",
    "10" = "SUBOPTIMAL", # ECOS_INACC_OFFSET
    "11" = "TIME_LIMIT", # Likely MI max iters reached
    "-1" = "TIME_LIMIT", # ECOS_MAXIT
    "-2" = "NUMERIC_ERROR"
  )
  exitflag <- as.character(result$retcodes["exitFlag"])
  status <- status_map[[exitflag]]
  if (is.null(status)) {
    status <- result$infostring
  }

  msg <- paste("Mixed Integer Linear Program Complete.")
  send_note_if(msg, show_progress, message)

  # Return results
  msg <- paste(
    "Status of IQR program:",
    status,
    "| Objective:",
    format(result$summary["pcost"], scientific = FALSE, digits = 10)
  )
  send_note_if(msg, !quietly, message)

  out$iqr <- list(
    obj = obj,
    G = G_ecos,
    h = h_ecos,
    A = A_ecos,
    b = b_ecos,
    bool_vars = bool_vars
  )
  out$params <- params
  out$result <- result
  out$status <- status

  if (status %in% c("OPTIMAL", "SUBOPTIMAL")) {
    answer <- result$x
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
    out$objval <- as.numeric(result$summary["pcost"])
  }

  return(out)
}
