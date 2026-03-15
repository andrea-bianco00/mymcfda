# =============================================================================
# kernel_h_sigma.R
# -----------------------------------------------------------------------------
# Scaled kernel function K_{h_sigma}(u) = (1/h_sigma) * K(u / h_sigma)
#
# This is the Block 2 analogue of kernel_h_mu.R. It is used in the local
# linear smoother for the estimation of varsigma^2(t), as described in
# Section 2.2 of Lin & Wang (2022):
#
#   varsigma2_hat(t) = argmin sum_i w_i sum_j K_{h_sigma}(T_ij - t)
#                      * [Z_ij - b0 - b1*(T_ij - t)]^2
#
# where Z_ij = {Y_ij - mu_hat(T_ij)}^2 are the squared residuals and
# K_{h_sigma}(T_ij - t) = (1/h_sigma) * K((T_ij - t) / h_sigma).
#
# Dependencies:
#   - kernel_fun()  defined in kernel.R
# =============================================================================


#' Scaled kernel function K_{h_sigma}(u) = (1/h_sigma) * K(u/h_sigma)
#'
#' @param u        Numeric vector. Raw argument, e.g. T_ij - t. Not rescaled.
#' @param h_sigma  Numeric scalar. Bandwidth h_sigma > 0.
#' @param kernel   Character string. Kernel type, passed to kernel_fun().
#'                 One of: "epanechnikov", "biweight", "triweight", "tricube",
#'                 "triangular", "cosine", "uniform", "gaussian".
#'                 No default — must be specified explicitly.
#'
#' @return Numeric vector of the same length as u with values K_{h_sigma}(u).
#'
#' @examples
#' u <- seq(-0.5, 0.5, by = 0.01)
#' plot(u, kernel_h_sigma(u, h_sigma = 0.3, kernel = "epanechnikov"), type = "l")
#' lines(u, kernel_h_sigma(u, h_sigma = 0.3, kernel = "biweight"), col = "red")

#' @export
kernel_h_sigma <- function(u, h_sigma, kernel) {
  
  # --- input checks -----------------------------------------------------------
  if (!is.numeric(u)) {
    stop("u must be a numeric vector.")
  }
  if (!is.numeric(h_sigma) || length(h_sigma) != 1 || h_sigma <= 0) {
    stop("h_sigma must be a single strictly positive number.")
  }
  
  # --- compute K_{h_sigma}(u) = (1/h_sigma) * K(u/h_sigma) ------------------
  (1/h_sigma) * kernel_fun(u/h_sigma, kernel = kernel)
}
