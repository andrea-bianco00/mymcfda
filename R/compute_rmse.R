# =============================================================================
# compute_rmse.R
# -----------------------------------------------------------------------------
# Computes RMSE for the noise variance estimator sigma^2_0, as defined
# in Section 5 of Lin & Wang (2022):
#
#   RMSE = sqrt( (1/N) * sum_{i=1}^{N} |sigma2_0_hat_i - sigma2_0|^2 )
#
# Dependencies: none (only base R)
# =============================================================================


#' Compute RMSE for sigma^2_0 estimates
#'
#' @param estimates     Numeric vector of length N. sigma^2_0_hat from each
#'                      simulation replicate.
#' @param sigma2_0_true Numeric scalar. The true noise variance sigma^2_0.
#' @param na_rm         Logical. If TRUE, removes NA/Inf values in estimates.
#'
#' @return Numeric scalar. The RMSE value.
#'
#' @examples
#' estimates <- c(0.11, 0.09, 0.13, 0.08, 0.12)
#' compute_rmse(estimates, sigma2_0_true = 0.10)

#' @export
compute_rmse <- function(estimates, sigma2_0_true, na_rm = FALSE) {
  
  # ---------------------------------------------------------------------------
  # input checks
  # ---------------------------------------------------------------------------
  
  if (!is.numeric(estimates) || length(estimates) == 0) {
    stop("'estimates' must be a non-empty numeric vector.")
  }
  
  if (!is.numeric(sigma2_0_true) || length(sigma2_0_true) != 1 ||
      is.na(sigma2_0_true) || !is.finite(sigma2_0_true)) {
    stop("'sigma2_0_true' must be a finite numeric scalar.")
  }
  
  if (!is.logical(na_rm) || length(na_rm) != 1) {
    stop("'na_rm' must be TRUE or FALSE.")
  }
  
  
  # ---------------------------------------------------------------------------
  # handle NA / Inf values
  # ---------------------------------------------------------------------------
  
  if (na_rm) {
    
    estimates <- estimates[is.finite(estimates)]
    
    if (length(estimates) == 0) {
      stop("After removing NA/Inf values, 'estimates' is empty.")
    }
    
  } else {
    
    if (any(is.na(estimates)) || any(!is.finite(estimates))) {
      stop("'estimates' contains NA or non-finite values. Set na_rm = TRUE if desired.")
    }
    
  }
  
  
  # ---------------------------------------------------------------------------
  # RMSE computation
  # ---------------------------------------------------------------------------
  
  sqrt(mean((estimates - sigma2_0_true)^2))
}
