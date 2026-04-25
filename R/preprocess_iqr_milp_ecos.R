### Meta -------------------------
###
### Title: Compute IQR estimator by preprocessing MILP using ECOS
###
### Description: The IQR estimator is computed by solving a MILP.
### This function pre-processes the data and fixes the sign of the residuals of
### outliers to speed up the procedure using the ECOS solver.
###
### Author: Omkar A. Katta (adapted for ECOS)
###

### preprocess_iqr_milp_ecos -------------------------
#' Compute IQR estimator by preprocessing MILP using ECOS
#'
#' @inheritParams preprocess_iqr_milp
#' @param ... Arguments passed to \code{iqr_milp_ecos}
#'
#' @export
preprocess_iqr_milp_ecos <- function(
  Y,
  D,
  X,
  Z,
  Phi = linear_projection(D, X, Z),
  tau,
  weights = NULL,
  M = NULL,
  TimeLimit = NULL,
  globalTimeLimit = Inf,
  prop_alpha_initial = 0.7,
  r = 1.25,
  show_iterations = FALSE,
  quietly = TRUE,
  show_progress = TRUE,
  ...
) {
  out <- list() # Initialize list of results to return

  # Start the clock
  clock_start <- Sys.time()

  # Get dimensions of data
  n <- length(Y)
  n_D <- nrow(D)
  n_X <- nrow(X)
  n_Z <- nrow(Z)
  p_D <- ncol(D)
  p_X <- ncol(X)
  p_Z <- ncol(Z)

  if (!is.null(weights)) {
    stopifnot(length(weights) == n)
    stopifnot(all(weights >= 0))
  }

  # If there are no endogeneous variables, return quantile regression results:
  if (p_D == 0) {
    msg <- paste("p_D is 0 -- running QR instead of IQR MILP...")
    send_note_if(msg, TRUE, warning)
    qr <- quantreg::rq(Y ~ X - 1, tau = tau, weights = weights)
    out$final_fit <- qr
    return(out)
  }

  # Determine preliminary residuals
  if (p_X == 0) {
    resid <- quantreg::rq(Y ~ D - 1, tau = tau, weights = weights)$residuals
  } else if (p_D == 0) {
    resid <- quantreg::rq(Y ~ X - 1, tau = tau, weights = weights)$residuals
  } else {
    resid <- quantreg::rq(Y ~ X + D - 1, tau = tau, weights = weights)$residuals
  }

  # Determine initial residual bounds
  alpha_initial <- stats::quantile(abs(resid), prop_alpha_initial)

  # Determine M before iqr_milp to avoid rerunning quantreg::rq
  if (is.null(M)) {
    max_qr <- max(abs(resid))
    M <- 2 * max_qr
  }

  # Start the while loop
  alphawidth <- alpha_initial
  status <- "TIME_LIMIT"
  num_fixed_vars_per_iteration <- c()
  time_limit_per_iteration <- c()
  time_elapsed_per_iteration <- c()
  final_objective_per_iteration <- c()

  counter <- 0
  obj <- 1 # initial value to enter loop
  while (status == "TIME_LIMIT" || obj > 1e-6) {
    counter <- counter + 1
    send_note_if(paste("Iteration", counter), show_iterations, message)
    while_start_time <- Sys.time()

    O_neg <- which(resid < -1 * alphawidth)
    O_pos <- which(resid > alphawidth)
    O <- c(O_neg, O_pos)
    send_note_if(paste("Alpha:", alphawidth), show_iterations, message)
    send_note_if(
      paste("Number of Fixed Dual Variables:", length(O)),
      show_iterations,
      message
    )
    num_fixed_vars_per_iteration <- c(num_fixed_vars_per_iteration, length(O))

    # Heuristic for time limit
    if (length(O) == 0) {
      TT <- Inf
    } else if (is.null(TimeLimit)) {
      num_free <- n - length(O)
      TT <- exp(num_free / 200 + p_D / 5 + num_free * p_D / 1000) * 4
    } else {
      TT <- TimeLimit
    }
    if (TT > globalTimeLimit) {
      TT <- globalTimeLimit
    }
    send_note_if(paste("TT:", TT), show_iterations, message)
    time_limit_per_iteration <- c(time_limit_per_iteration, TT)

    # IQR using ECOS
    fit <- iqr_milp_ecos(
      Y = Y,
      X = X,
      D = D,
      Z = Z,
      Phi = Phi,
      tau = tau,
      weights = weights,
      O_neg = O_neg,
      O_pos = O_pos,
      TimeLimit = TT,
      M = M,
      quietly = quietly,
      show_progress = show_progress,
      ...
    )

    final_objective_per_iteration <- c(
      final_objective_per_iteration,
      ifelse(is.null(fit$objval), "NULL", as.character(round(fit$objval, 3)))
    )
    if (is.null(fit$objval)) {
      obj <- 0.5
    } else {
      obj <- fit$objval
    }
    status <- fit$status
    alphawidth <- alphawidth * r
    send_note_if(
      paste("Iteration", counter, "complete"),
      show_iterations,
      message
    )
    if (TT == Inf & obj > 1e-6) {
      warning("Nonzero Coefficients on Instruments")
      break
    }
    current <- Sys.time()
    while_elapsed <- difftime(current, while_start_time, units = "secs")
    time_elapsed_per_iteration <- c(time_elapsed_per_iteration, while_elapsed)

    elapsed_time <- difftime(current, clock_start, units = "secs")
    if (as.numeric(elapsed_time) > globalTimeLimit) {
      warning(paste("Global Time Limit of", globalTimeLimit, "reached."))
      break
    }
  }

  # Stop the clock
  clock_end <- Sys.time()
  elapsed_time <- difftime(clock_end, clock_start, units = "mins")

  # Return results
  out$status <- fit$status
  out$final_fit <- fit
  out$time <- elapsed_time # mins
  out$O_neg <- O_neg
  out$O_pos <- O_pos
  out$iterations <- counter
  out$num_fixed_vars_per_iteration <- num_fixed_vars_per_iteration
  out$time_limit_per_iteration <- time_limit_per_iteration
  out$time_elapsed_per_iteration <- time_elapsed_per_iteration
  out$final_objective_per_iteration <- final_objective_per_iteration

  return(out)
}
