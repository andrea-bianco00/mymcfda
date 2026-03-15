# =============================================================================
# varsigma2_hat.R
# -----------------------------------------------------------------------------
# Computes the ridged local linear estimator of varsigma^2(t) on a grid of
# points t_1, ..., t_k.
#
# This is the Block 2 analogue of mu_hat.R. It is a wrapper around
# compute_varsigma2_hat_t(). Since compute_varsigma2_hat_t() now accepts a
# vector of evaluation points, the entire grid is computed in a single call —
# no loop or sapply needed.
#
# If t_grid is not provided, it is built automatically from the data as a
# uniform grid over [min(T_ij), max(T_ij)].
#
# A non-negativity correction is optionally applied at the wrapper level:
# since varsigma^2(t) = sigma_X^2(t) + sigma_0^2 >= 0 by definition, any
# negative estimate is treated as a numerical artifact and truncated to 0.
#
# Dependencies:
#   - compute_varsigma2_hat_t()  defined in varsigma2_hat_t.R
#   - compute_Sr_sigma()         defined in S_r_sigma.R
#   - compute_Rr_sigma()         defined in R_r_sigma.R
#   - compute_weights()          defined in weights.R
#   - kernel_h_sigma()           defined in kernel_h_sigma.R
#   - kernel_fun()               defined in kernel.R
# =============================================================================


#' Ridged local linear estimator of varsigma^2(t) on a grid of points
#'
#' @param T_list                List of length n. T_list[[i]] is the vector of
#'                              observation times for subject i.
#' @param Z_list                List of length n. Z_list[[i]] is the vector of
#'                              squared residuals Z_ij = {Y_ij - mu_hat(T_ij)}^2
#'                              for subject i.
#' @param h_sigma               Numeric scalar. Bandwidth h_sigma > 0.
#' @param kernel                Character string. Kernel type, passed to
#'                              compute_varsigma2_hat_t().
#' @param scheme                Character string. Weighting scheme.
#'                              One of "OBS", "SUBJ", "OPTIMAL".
#' @param t_grid                Numeric vector. Grid of points at which
#'                              varsigma^2(t) is estimated. If NULL (default),
#'                              built automatically as n_t equally spaced points
#'                              over [min(T_ij), max(T_ij)].
#' @param n_t                   Positive integer. Number of grid points when
#'                              t_grid is built automatically. Default is 100.
#' @param truncate_nonnegative  Logical. If TRUE (default), truncate negative
#'                              estimates to 0.
#'
#' @return A list with:
#'   $t_grid:        Numeric vector. The grid of evaluation points.
#'   $varsigma2_hat: Numeric vector. Estimated varsigma2_hat(t) for each t in t_grid.

compute_varsigma2_hat <- function(T_list, Z_list, h_sigma, kernel, scheme,
                                  t_grid = NULL, n_t = 100,
                                  truncate_nonnegative = TRUE) {
  
  # --- input checks -----------------------------------------------------------
  if (!is.list(T_list) || !is.list(Z_list))
    stop("T_list and Z_list must be lists.")
  if (length(T_list) != length(Z_list))
    stop("T_list and Z_list must have the same length n.")
  
  # --- build t_grid if not provided -------------------------------------------
  if (is.null(t_grid)) {
    t_min  <- min(sapply(T_list, min))
    t_max  <- max(sapply(T_list, max))
    t_grid <- seq(t_min, t_max, length.out = n_t)
    cat(sprintf("t_grid built automatically: %d points in [%.4f, %.4f]\n",
                n_t, t_min, t_max))
  } else {
    if (!is.numeric(t_grid) || length(t_grid) < 1)
      stop("t_grid must be a non-empty numeric vector.")
  }
  
  # --- evaluate varsigma2_hat on the entire grid in one call ------------------
  varsigma2_hat <- compute_varsigma2_hat_t(
    t_vec   = t_grid,
    T_list  = T_list,
    Z_list  = Z_list,
    h_sigma = h_sigma,
    kernel  = kernel,
    scheme  = scheme
  )
  
  # --- non-negativity correction ----------------------------------------------
  if (truncate_nonnegative) {
    varsigma2_hat <- pmax(varsigma2_hat, 0)
  }
  
  return(list(
    t_grid        = t_grid,
    varsigma2_hat = varsigma2_hat
  ))
}
