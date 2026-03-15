# =============================================================================
# compute_rmise_covariance.R
# -----------------------------------------------------------------------------
# Computes the RMISE for the covariance function C(s,t):
#
#   RMISE = sqrt( (1/N) * sum_{i=1}^N \int_T \int_T
#                    (C_hat_i(s,t) - C_true(s,t))^2 ds dt )
#
# Numerical integration is performed by the 2D trapezoidal rule on the
# rectangular grid s_grid x t_grid.
# =============================================================================


# -----------------------------------------------------------------------------
# Helper: trapezoidal weights on a 1D grid
# -----------------------------------------------------------------------------
trapz_weights_1d <- function(x) {
  
  if (!is.numeric(x) || length(x) < 2) {
    stop("'x' must be a numeric vector with at least two points.")
  }
  
  if (any(is.na(x))) {
    stop("'x' must not contain NA values.")
  }
  
  if (any(diff(x) <= 0)) {
    stop("'x' must be strictly increasing.")
  }
  
  n <- length(x)
  dx <- diff(x)
  w  <- numeric(n)
  
  w[1] <- dx[1] / 2
  w[n] <- dx[n - 1] / 2
  
  if (n > 2) {
    for (k in 2:(n - 1)) {
      w[k] <- (dx[k - 1] + dx[k]) / 2
    }
  }
  
  w
}


# -----------------------------------------------------------------------------
# Helper: 2D trapezoidal integral
# -----------------------------------------------------------------------------
trapz_2d <- function(x, y, z_mat) {
  
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("'x' and 'y' must be numeric.")
  }
  
  if (!is.matrix(z_mat)) {
    stop("'z_mat' must be a matrix.")
  }
  
  if (any(is.na(x)) || any(is.na(y)) || any(is.na(z_mat))) {
    stop("'x', 'y', and 'z_mat' must not contain NA values.")
  }
  
  if (any(diff(x) <= 0)) {
    stop("'x' must be strictly increasing.")
  }
  
  if (any(diff(y) <= 0)) {
    stop("'y' must be strictly increasing.")
  }
  
  if (!all(dim(z_mat) == c(length(x), length(y)))) {
    stop("'z_mat' must have dimensions length(x) x length(y).")
  }
  
  wx <- trapz_weights_1d(x)
  wy <- trapz_weights_1d(y)
  
  sum(outer(wx, wy) * z_mat)
}


# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
#' Compute RMISE for the covariance function C(s,t)
#'
#' @param estimates_list  List of length N. Each element is either:
#'                        - a matrix giving C_hat(s,t) on s_grid x t_grid, or
#'                        - a list containing component $covariance_hat.
#' @param s_grid          Numeric vector. Grid for the first argument s.
#' @param t_grid          Numeric vector. Grid for the second argument t.
#' @param C_true_mat      Optional matrix. True covariance on s_grid x t_grid.
#' @param C_true_fun      Optional function of two arguments (s, t), returning
#'                        the true covariance matrix on the grid.
#'
#' Exactly one of C_true_mat or C_true_fun must be provided.
#'
#' @return Numeric scalar. The RMISE value.
#'
#' @examples
#' s_grid <- seq(0, 1, length.out = 50)
#' t_grid <- seq(0, 1, length.out = 50)
#'
#' C_true_fun <- function(s, t) {
#'   outer(s, t, function(x, y) exp(-abs(x - y)))
#' }
#'
#' C1 <- C_true_fun(s_grid, t_grid) + matrix(rnorm(50 * 50, 0, 0.01), 50, 50)
#' C2 <- C_true_fun(s_grid, t_grid) + matrix(rnorm(50 * 50, 0, 0.02), 50, 50)
#'
#' compute_rmise_covariance(
#'   estimates_list = list(C1, C2),
#'   s_grid         = s_grid,
#'   t_grid         = t_grid,
#'   C_true_fun     = C_true_fun
#' )

compute_rmise_covariance <- function(estimates_list,
                                     s_grid,
                                     t_grid,
                                     C_true_mat = NULL,
                                     C_true_fun = NULL) {
  
  # ---------------------------------------------------------------------------
  # 1. input checks
  # ---------------------------------------------------------------------------
  if (!is.list(estimates_list) || length(estimates_list) == 0) {
    stop("'estimates_list' must be a non-empty list.")
  }
  
  if (!is.numeric(s_grid) || length(s_grid) < 2) {
    stop("'s_grid' must be a numeric vector with at least two points.")
  }
  
  if (!is.numeric(t_grid) || length(t_grid) < 2) {
    stop("'t_grid' must be a numeric vector with at least two points.")
  }
  
  if (any(is.na(s_grid)) || any(is.na(t_grid))) {
    stop("'s_grid' and 't_grid' must not contain NA values.")
  }
  
  if (any(diff(s_grid) <= 0)) {
    stop("'s_grid' must be strictly increasing.")
  }
  
  if (any(diff(t_grid) <= 0)) {
    stop("'t_grid' must be strictly increasing.")
  }
  
  if (!is.null(C_true_mat) && !is.null(C_true_fun)) {
    stop("Provide only one of 'C_true_mat' or 'C_true_fun', not both.")
  }
  
  if (is.null(C_true_mat) && is.null(C_true_fun)) {
    stop("You must provide either 'C_true_mat' or 'C_true_fun'.")
  }
  
  # ---------------------------------------------------------------------------
  # 2. evaluate true covariance on grid
  # ---------------------------------------------------------------------------
  if (!is.null(C_true_fun)) {
    
    if (!is.function(C_true_fun)) {
      stop("'C_true_fun' must be a function.")
    }
    
    C_true_vals <- C_true_fun(s_grid, t_grid)
    
  } else {
    
    C_true_vals <- C_true_mat
  }
  
  if (!is.matrix(C_true_vals)) {
    stop("True covariance must be a matrix.")
  }
  
  if (!all(dim(C_true_vals) == c(length(s_grid), length(t_grid)))) {
    stop("True covariance matrix must have dimensions length(s_grid) x length(t_grid).")
  }
  
  if (any(is.na(C_true_vals)) || any(!is.finite(C_true_vals))) {
    stop("True covariance matrix must be finite and must not contain NA values.")
  }
  
  # ---------------------------------------------------------------------------
  # 3. compute integrated squared error for each replicate
  # ---------------------------------------------------------------------------
  ise_vec <- numeric(length(estimates_list))
  
  for (i in seq_along(estimates_list)) {
    
    est_i <- estimates_list[[i]]
    
    # allow either a raw matrix or a list with component $covariance_hat
    if (is.list(est_i)) {
      if (!("covariance_hat" %in% names(est_i))) {
        stop(sprintf(
          "estimates_list[[%d]] is a list but does not contain 'covariance_hat'.",
          i
        ))
      }
      est_i <- est_i$covariance_hat
    }
    
    if (!is.matrix(est_i)) {
      stop(sprintf("Replicate %d: estimated covariance must be a matrix.", i))
    }
    
    if (!all(dim(est_i) == c(length(s_grid), length(t_grid)))) {
      stop(sprintf(
        "Replicate %d: estimated covariance must have dimensions length(s_grid) x length(t_grid).",
        i
      ))
    }
    
    if (any(is.na(est_i)) || any(!is.finite(est_i))) {
      stop(sprintf(
        "Replicate %d: estimated covariance contains NA or non-finite values.",
        i
      ))
    }
    
    sq_error_i <- (est_i - C_true_vals)^2
    
    ise_vec[i] <- trapz_2d(s_grid, t_grid, sq_error_i)
  }
  
  # ---------------------------------------------------------------------------
  # 4. RMISE
  # ---------------------------------------------------------------------------
  sqrt(mean(ise_vec))
}
