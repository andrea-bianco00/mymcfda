# =============================================================================
# mu_hat.R
# -----------------------------------------------------------------------------
# Computes the ridged local linear estimator of the mean function mu(t)
# on a grid of points t_1, ..., t_k.
#
# This function is a wrapper around compute_mu_hat_t(). Since compute_mu_hat_t()
# now accepts a vector of evaluation points, the entire grid is computed in a
# single call — no loop or sapply needed.
#
# If t_grid is not provided, it is built automatically from the data as a
# uniform grid over [min(T_ij), max(T_ij)].
#
# Dependencies:
#   - compute_mu_hat_t()  defined in mu_hat_t.R
#   - compute_Sr()        defined in S_r.R
#   - compute_Rr()        defined in R_r.R
#   - compute_weights()   defined in weights.R
#   - kernel_h_mu()       defined in kernel_h_mu.R
#   - kernel_fun()        defined in kernel.R
# =============================================================================


#' Ridged local linear estimator of mu(t) on a grid of points
#'
#' @param T_list  List of length n. T_list((i)) is the vector of observation
#'                times for subject i.
#' @param Y_list  List of length n. Y_list((i)) is the vector of observed
#'                values Y_ij for subject i.
#' @param h_mu    Numeric scalar. Bandwidth h_mu > 0.
#' @param kernel  Character string. Kernel type, passed to compute_mu_hat_t().
#' @param scheme  Character string. Weighting scheme. One of "OBS", "SUBJ", "OPTIMAL".
#' @param t_grid  Numeric vector. Grid of points at which mu(t) is estimated.
#'                If NULL (default), built automatically as n_t equally spaced
#'                points over (min(T_ij), max(T_ij)).
#' @param n_t     Positive integer. Number of grid points when t_grid is built
#'                automatically. Default is 100.
#'
#' @return A list with:
#'   $t_grid:  Numeric vector. The grid of evaluation points.
#'   $mu_hat:  Numeric vector. Estimated mu_hat(t) for each t in t_grid.

#' @export
compute_mu_hat <- function(T_list, Y_list, h_mu, kernel, scheme,
                           t_grid = NULL, n_t = 100) {
  
  # --- input checks -----------------------------------------------------------
  if (!is.list(T_list) || !is.list(Y_list))
    stop("T_list and Y_list must be lists.")
  if (length(T_list) != length(Y_list))
    stop("T_list and Y_list must have the same length n.")
  
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
  
  # --- evaluate mu_hat on the entire grid in one call -------------------------
  mu_hat <- compute_mu_hat_t(t_vec  = t_grid,
                             T_list = T_list,
                             Y_list = Y_list,
                             h_mu   = h_mu,
                             kernel = kernel,
                             scheme = scheme)
  
  return(list(t_grid = t_grid, mu_hat = mu_hat))
}
