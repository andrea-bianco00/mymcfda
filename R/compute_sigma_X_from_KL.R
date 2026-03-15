# =============================================================================
# compute_sigma_X_from_KL.R
# -----------------------------------------------------------------------------
# Computes the variance function sigma_X^2(t) as a CONSEQUENCE of the
# Karhunen-Loève specification:
#
#   sigma_X^2(t) = C(t, t) = sum_{k=1}^{K} lambda_k * phi_k(t)^2
#
# This function does NOT require the user to specify sigma_X^2(t) — it
# derives it automatically from the eigenfunctions and eigenvalues.
# In the KL strada, sigma_X^2(t) is not an input but a consequence.
#
# PURPOSE:
#   - Diagnostics: visualise the variance of the generated process
#   - Comparison: compare with the estimated sigma_X^2(t) from the
#     estimation pipeline to assess estimator performance
#   - Verification: check that the KL specification produces a
#     reasonable variance structure
#
# INPUT:
#   t             : numeric vector of time points where to evaluate
#                   sigma_X^2(t). Typically a fine grid over [0, 1].
#   eigenfn_list  : list of K vectorised functions phi_k(t).
#                   Same as passed to generate_KL_process().
#   eigenval_vec  : numeric vector of K positive eigenvalues lambda_k.
#                   Same as passed to generate_KL_process().
#
# OUTPUT:
#   Numeric vector of the same length as t, containing
#   sigma_X^2(t_1), ..., sigma_X^2(t_G).
#
# EXAMPLES:
#   # Cov II from Lin & Wang (2022)
#   K <- 50
#   eigenfn_list <- lapply(1:K, function(k) {
#     function(t) sqrt(2) * sin(2 * k * pi * t)
#   })
#   eigenval_vec <- 2 * (1:K)^(-2)
#
#   t_grid <- seq(0, 1, length.out = 200)
#   sigma2 <- compute_sigma_X_from_KL(t_grid, eigenfn_list, eigenval_vec)
#   plot(t_grid, sigma2, type = "l", main = "Variance function from KL")
#
# Dependencies: base R only.
# =============================================================================


#' Compute variance function from KL specification
#'
#' @param t             Numeric vector of evaluation points.
#' @param eigenfn_list  List of K vectorised functions phi_k(t).
#' @param eigenval_vec  Numeric vector of K positive eigenvalues.
#'
#' @return Numeric vector of same length as t with sigma_X^2(t) values.

compute_sigma_X_from_KL <- function(t, eigenfn_list, eigenval_vec) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(t) || length(t) < 1)
    stop("t must be a non-empty numeric vector.")
  
  if (!is.list(eigenfn_list) || length(eigenfn_list) < 1)
    stop("eigenfn_list must be a non-empty list of functions.")
  
  K <- length(eigenfn_list)
  
  if (!is.numeric(eigenval_vec) || length(eigenval_vec) != K)
    stop("eigenval_vec must be a numeric vector of length K = ", K, ".")
  
  if (any(eigenval_vec <= 0))
    stop("All eigenvalues must be strictly positive.")
  
  # ---------------------------------------------------------------------------
  # COMPUTE sigma_X^2(t) = sum_k lambda_k * phi_k(t)^2
  # ---------------------------------------------------------------------------
  sigma2 <- rep(0, length(t))
  
  for (k in seq_len(K)) {
    
    if (!is.function(eigenfn_list[[k]]))
      stop("eigenfn_list[[", k, "]] is not a function.")
    
    phi_k_vals <- eigenfn_list[[k]](t)
    
    if (!is.numeric(phi_k_vals) || length(phi_k_vals) != length(t))
      stop("eigenfn_list[[", k, "]] must be vectorised: it must return a ",
           "numeric vector of the same length as t.")
    
    sigma2 <- sigma2 + eigenval_vec[k] * phi_k_vals^2
  }
  
  return(sigma2)
}
