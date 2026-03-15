# =============================================================================
# compute_rmise_variance.R
# -----------------------------------------------------------------------------
# Computes the RMISE for the variance function sigma_X^2(t):
#
#   RMISE = sqrt( (1/N) * sum_{i=1}^N integral_T
#                    (sigmaX2_hat_i(t) - sigmaX2_true(t))^2 dt )
#
# Numerical integration is performed by the trapezoidal rule on t_grid.
# =============================================================================


# -----------------------------------------------------------------------------
# Helper: 1D trapezoidal integral
# -----------------------------------------------------------------------------
trapz_1d <- function(x, y) {
  
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("'x' and 'y' must be numeric.")
  }
  
  if (length(x) != length(y)) {
    stop("'x' and 'y' must have the same length.")
  }
  
  if (length(x) < 2) {
    stop("'x' must contain at least two points.")
  }
  
  if (any(is.na(x)) || any(is.na(y))) {
    stop("'x' and 'y' must not contain NA values.")
  }
  
  if (any(diff(x) <= 0)) {
    stop("'x' must be strictly increasing.")
  }
  
  sum((y[-1] + y[-length(y)]) / 2 * diff(x))
}


# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
#' Compute RMISE for sigma_X^2(t)
#'
#' @param estimates_list     List of length N. Each element is either:
#'                           - a numeric vector giving sigmaX2_hat on t_grid, or
#'                           - a list containing component $sigmaX2_hat.
#' @param t_grid             Numeric vector. Common evaluation grid.
#' @param sigmaX2_true       Optional numeric vector. True sigma_X^2(t) on t_grid.
#' @param sigmaX2_true_fun   Optional function. True sigma_X^2(t), evaluated on t_grid.
#'
#' Exactly one of sigmaX2_true or sigmaX2_true_fun must be provided.
#'
#' @return Numeric scalar. The RMISE value.
#'
#' @examples
#' t_grid <- seq(0, 1, length.out = 100)
#' sigmaX2_true_fun <- function(t) 0.2 + 0.1 * sin(2*pi*t)
#'
#' est1 <- sigmaX2_true_fun(t_grid) + rnorm(100, 0, 0.01)
#' est2 <- sigmaX2_true_fun(t_grid) + rnorm(100, 0, 0.02)
#'
#' compute_rmise_variance(
#'   estimates_list   = list(est1, est2),
#'   t_grid           = t_grid,
#'   sigmaX2_true_fun = sigmaX2_true_fun
#' )

compute_rmise_variance <- function(estimates_list,
                                   t_grid,
                                   sigmaX2_true = NULL,
                                   sigmaX2_true_fun = NULL) {
  
  # ---------------------------------------------------------------------------
  # 1. input checks
  # ---------------------------------------------------------------------------
  if (!is.list(estimates_list) || length(estimates_list) == 0) {
    stop("'estimates_list' must be a non-empty list.")
  }
  
  if (!is.numeric(t_grid) || length(t_grid) < 2) {
    stop("'t_grid' must be a numeric vector with at least two points.")
  }
  
  if (any(is.na(t_grid))) {
    stop("'t_grid' must not contain NA values.")
  }
  
  if (any(diff(t_grid) <= 0)) {
    stop("'t_grid' must be strictly increasing.")
  }
  
  if (!is.null(sigmaX2_true) && !is.null(sigmaX2_true_fun)) {
    stop("Provide only one of 'sigmaX2_true' or 'sigmaX2_true_fun', not both.")
  }
  
  if (is.null(sigmaX2_true) && is.null(sigmaX2_true_fun)) {
    stop("You must provide either 'sigmaX2_true' or 'sigmaX2_true_fun'.")
  }
  
  # ---------------------------------------------------------------------------
  # 2. evaluate true sigmaX2 on grid
  # ---------------------------------------------------------------------------
  if (!is.null(sigmaX2_true_fun)) {
    
    if (!is.function(sigmaX2_true_fun)) {
      stop("'sigmaX2_true_fun' must be a function.")
    }
    
    sigmaX2_true_vals <- sigmaX2_true_fun(t_grid)
    
  } else {
    
    sigmaX2_true_vals <- sigmaX2_true
  }
  
  if (!is.numeric(sigmaX2_true_vals) || length(sigmaX2_true_vals) != length(t_grid)) {
    stop("True sigmaX2 values must be a numeric vector of the same length as 't_grid'.")
  }
  
  if (any(is.na(sigmaX2_true_vals)) || any(!is.finite(sigmaX2_true_vals))) {
    stop("True sigmaX2 values must be finite and must not contain NA.")
  }
  
  # ---------------------------------------------------------------------------
  # 3. compute integrated squared error for each replicate
  # ---------------------------------------------------------------------------
  ise_vec <- numeric(length(estimates_list))
  
  for (i in seq_along(estimates_list)) {
    
    est_i <- estimates_list[[i]]
    
    # allow either a raw vector or a list with component $sigmaX2_hat
    if (is.list(est_i)) {
      if (!("sigmaX2_hat" %in% names(est_i))) {
        stop(sprintf(
          "estimates_list[[%d]] is a list but does not contain 'sigmaX2_hat'.",
          i
        ))
      }
      est_i <- est_i$sigmaX2_hat
    }
    
    if (!is.numeric(est_i) || length(est_i) != length(t_grid)) {
      stop(sprintf(
        "Replicate %d: estimated sigmaX2 must be numeric and have length equal to length(t_grid).",
        i
      ))
    }
    
    if (any(is.na(est_i)) || any(!is.finite(est_i))) {
      stop(sprintf(
        "Replicate %d: estimated sigmaX2 contains NA or non-finite values.",
        i
      ))
    }
    
    sq_error_i <- (est_i - sigmaX2_true_vals)^2
    
    ise_vec[i] <- trapz_1d(t_grid, sq_error_i)
  }
  
  # ---------------------------------------------------------------------------
  # 4. RMISE
  # ---------------------------------------------------------------------------
  sqrt(mean(ise_vec))
}
