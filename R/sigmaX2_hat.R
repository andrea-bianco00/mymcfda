# =============================================================================
# sigmaX2_hat.R
# -----------------------------------------------------------------------------
# Computes the final estimator of sigma_X^2(t):
#
#   sigma_X^2_hat(t) = varsigma^2_hat(t) - sigma_0^2_hat
#
# Since sigma_X^2(t) >= 0 by definition, negative values are truncated to zero.
#
# This function returns both:
#   - the raw estimator before truncation
#   - the final nonnegative estimator after truncation
#
# Dependencies:
#   - res_varsigma must be the output of compute_varsigma2_hat()
#   - sigma0_hat must be a scalar estimate of sigma_0^2
# =============================================================================


#' Compute the final estimator of sigma_X^2(t)
#'
#' @param res_varsigma List. Output of compute_varsigma2_hat(), containing:
#'   - t_grid
#'   - varsigma2_hat
#' @param sigma0_hat Numeric scalar. Estimated noise variance sigma_0^2.
#' @param truncate_nonnegative Logical. If TRUE (default), truncate negative
#'   values of sigmaX2_raw to zero.
#' @param warn Logical. If TRUE (default), issue a warning when truncation occurs.
#'
#' @return A list with:
#'   - t_grid:        grid of evaluation points
#'   - sigmaX2_raw:   raw estimator varsigma2_hat - sigma0_hat
#'   - sigmaX2_hat:   final nonnegative estimator
#'   - n_truncated:   number of grid points truncated to zero
#'
#' @examples
#' # res_varsigma <- compute_varsigma2_hat(...)
#' # sigma0_res   <- compute_sigma0(...)
#' # res_sigmaX2  <- compute_sigmaX2_hat(res_varsigma, sigma0_res$sigma0_hat)

compute_sigmaX2_hat <- function(res_varsigma,
                                sigma0_hat,
                                truncate_nonnegative = TRUE,
                                warn = TRUE) {
  
  # --- input checks -----------------------------------------------------------
  if (!is.list(res_varsigma)) {
    stop("res_varsigma must be a list.")
  }
  
  required_names <- c("t_grid", "varsigma2_hat")
  if (!all(required_names %in% names(res_varsigma))) {
    stop("res_varsigma must contain 't_grid' and 'varsigma2_hat'.")
  }
  
  if (!is.numeric(sigma0_hat) || length(sigma0_hat) != 1 || is.na(sigma0_hat)) {
    stop("sigma0_hat must be a non-missing numeric scalar.")
  }
  
  t_grid <- res_varsigma$t_grid
  varsigma2_hat <- res_varsigma$varsigma2_hat
  
  if (!is.numeric(t_grid) || !is.numeric(varsigma2_hat)) {
    stop("'t_grid' and 'varsigma2_hat' must be numeric.")
  }
  
  if (length(t_grid) != length(varsigma2_hat)) {
    stop("'t_grid' and 'varsigma2_hat' must have the same length.")
  }
  
  # --- raw estimator ----------------------------------------------------------
  sigmaX2_raw <- varsigma2_hat - sigma0_hat
  
  # --- final estimator --------------------------------------------------------
  if (truncate_nonnegative) {
    sigmaX2_hat <- pmax(sigmaX2_raw, 0)
  } else {
    sigmaX2_hat <- sigmaX2_raw
  }
  
  # --- truncation diagnostics -------------------------------------------------
  n_truncated <- sum(sigmaX2_raw < 0)
  
  if (warn && truncate_nonnegative && n_truncated > 0) {
    warning(sprintf(
      "%d grid point(s) had varsigma2_hat(t) < sigma0_hat and were truncated to 0. This may indicate that sigma0_hat is overestimated or h_sigma is too large.",
      n_truncated
    ))
  }
  
  # --- output -----------------------------------------------------------------
  return(list(
    t_grid       = t_grid,
    sigmaX2_raw  = sigmaX2_raw,
    sigmaX2_hat  = sigmaX2_hat,
    n_truncated  = n_truncated
  ))
}
