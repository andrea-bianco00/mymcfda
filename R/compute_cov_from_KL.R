# =============================================================================
# compute_cov_from_KL.R
# -----------------------------------------------------------------------------
# Computes the covariance function C(s, t) as a CONSEQUENCE of the
# Karhunen-Loève specification:
#
#   C(s, t) = sum_{k=1}^{K} lambda_k * phi_k(s) * phi_k(t)
#
# This function does NOT require the user to specify C(s, t) — it
# derives it automatically from the eigenfunctions and eigenvalues.
# In the KL strada, C(s, t) is not an input but a consequence.
#
# PURPOSE:
#   - Diagnostics: visualise the true covariance surface
#   - Comparison: compare with the estimated C_hat(s, t) from the
#     estimation pipeline to assess estimator performance (e.g. RMISE)
#   - Verification: check that the KL specification produces a
#     valid covariance structure
#
# NOTE: when s = t, C(t, t) = sigma_X^2(t) = sum_k lambda_k * phi_k(t)^2,
#       which is the same quantity computed by compute_sigma_X_from_KL().
#
# INPUT:
#   s             : numeric vector of length G_s. First set of time points.
#   t             : numeric vector of length G_t. Second set of time points.
#   eigenfn_list  : list of K vectorised functions phi_k(t).
#                   Same as passed to generate_KL_process().
#   eigenval_vec  : numeric vector of K positive eigenvalues lambda_k.
#                   Same as passed to generate_KL_process().
#
# OUTPUT:
#   Numeric matrix of dimension G_s x G_t. Entry (i, j) contains
#   C(s_i, t_j) = sum_k lambda_k * phi_k(s_i) * phi_k(t_j).
#
# EXAMPLES:
#   # Cov II from Lin & Wang (2022)
#   K <- 50
#   eigenfn_list <- lapply(1:K, function(k) {
#     function(t) sqrt(2) * sin(2 * k * pi * t)
#   })
#   eigenval_vec <- 2 * (1:K)^(-2)
#
#   t_grid <- seq(0, 1, length.out = 100)
#   C_mat <- compute_cov_from_KL(t_grid, t_grid, eigenfn_list, eigenval_vec)
#   image(t_grid, t_grid, C_mat, main = "True covariance surface")
#
# Dependencies: base R only.
# =============================================================================


#' Compute covariance function from KL specification
#'
#' @param s             Numeric vector of first set of evaluation points.
#' @param t             Numeric vector of second set of evaluation points.
#' @param eigenfn_list  List of K vectorised functions phi_k(t).
#' @param eigenval_vec  Numeric vector of K positive eigenvalues.
#'
#' @return Numeric matrix of dimension length(s) x length(t) with
#'         C(s_i, t_j) values.

#' @export
compute_cov_from_KL <- function(s, t, eigenfn_list, eigenval_vec) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(s) || length(s) < 1)
    stop("s must be a non-empty numeric vector.")
  
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
  # COMPUTE C(s, t) = sum_k lambda_k * phi_k(s) * phi_k(t)
  #
  # Strategy: build two matrices
  #   Phi_s : G_s x K, with Phi_s[i, k] = sqrt(lambda_k) * phi_k(s_i)
  #   Phi_t : G_t x K, with Phi_t[j, k] = sqrt(lambda_k) * phi_k(t_j)
  # Then C = Phi_s %*% t(Phi_t), which is G_s x G_t.
  # This is efficient and avoids a triple loop.
  # ---------------------------------------------------------------------------
  G_s <- length(s)
  G_t <- length(t)
  sqrt_eigenvals <- sqrt(eigenval_vec)
  
  Phi_s <- matrix(NA_real_, nrow = G_s, ncol = K)
  Phi_t <- matrix(NA_real_, nrow = G_t, ncol = K)
  
  for (k in seq_len(K)) {
    
    if (!is.function(eigenfn_list[[k]]))
      stop("eigenfn_list[[", k, "]] is not a function.")
    
    phi_k_s <- eigenfn_list[[k]](s)
    phi_k_t <- eigenfn_list[[k]](t)
    
    if (!is.numeric(phi_k_s) || length(phi_k_s) != G_s)
      stop("eigenfn_list[[", k, "]] must be vectorised. ",
           "Expected length ", G_s, " for s, got ", length(phi_k_s), ".")
    
    if (!is.numeric(phi_k_t) || length(phi_k_t) != G_t)
      stop("eigenfn_list[[", k, "]] must be vectorised. ",
           "Expected length ", G_t, " for t, got ", length(phi_k_t), ".")
    
    Phi_s[, k] <- sqrt_eigenvals[k] * phi_k_s
    Phi_t[, k] <- sqrt_eigenvals[k] * phi_k_t
  }
  
  # C = Phi_s %*% t(Phi_t)  gives (G_s x K) %*% (K x G_t) = (G_s x G_t)
  C_mat <- Phi_s %*% t(Phi_t)
  
  return(C_mat)
}
