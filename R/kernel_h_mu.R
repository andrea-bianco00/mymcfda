# =============================================================================
# kernel_h_mu.R
# -----------------------------------------------------------------------------
# Scaled kernel function K_{h_mu}(u) = (1/h_mu) * K(u / h_mu)
#
# This function builds directly on kernel_fun() defined in kernel.R, which
# must be sourced before using this file.
#
# The scaled kernel K_{h_mu} is the version actually used in the local linear
# smoother for mean estimation (equation 2, Lin & Wang 2022):
#
#   mu_hat(t) = argmin sum_i w_i sum_j K_{h_mu}(T_ij - t) * {Y_ij - b0 - b1*(T_ij - t)}^2
#
# where K_{h_mu}(T_ij - t) = (1/h_mu) * K((T_ij - t) / h_mu)
# =============================================================================


#' Scaled kernel function K_{h_mu}(u) = (1/h_mu) * K(u/h_mu)
#'
#' @param u       Numeric vector. Raw argument, e.g. T_ij - t. Not rescaled.
#' @param h_mu    Numeric scalar. Bandwidth h_mu > 0.
#' @param kernel  Character string. Kernel type, passed to kernel_fun().
#'                One of: "epanechnikov", "biweight", "triweight", "tricube",
#'                "triangular", "cosine", "uniform", "gaussian".
#'                No default — must be specified explicitly.
#'
#' @return Numeric vector of the same length as u with values K_{h_mu}(u).
#'
#' @examples
#' u <- seq(-0.5, 0.5, by = 0.01)
#' plot(u, kernel_h_mu(u, h_mu = 0.3, kernel = "epanechnikov"), type = "l")
#' lines(u, kernel_h_mu(u, h_mu = 0.3, kernel = "biweight"), col = "red")

kernel_h_mu <- function(u, h_mu, kernel) {
  
  # --- input checks -----------------------------------------------------------
  if (!is.numeric(u)) {
    stop("u must be a numeric vector.")
  }
  if (!is.numeric(h_mu) || length(h_mu) != 1 || h_mu <= 0) {
    stop("h_mu must be a single strictly positive number.")
  }
  
  # --- compute K_{h_mu}(u) = (1/h_mu) * K(u/h_mu) ----------------------------
  (1/h_mu) * kernel_fun(u/h_mu, kernel = kernel)
}
