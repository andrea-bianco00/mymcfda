# =============================================================================
# R_r_sigma.R
# -----------------------------------------------------------------------------
# Computes the quantity R_r(t) used in the local linear estimator of
# varsigma^2(t) in Block 2.
#
# This is the Block 2 analogue of R_r.R. The formula is identical, but:
#   - uses h_sigma instead of h_mu
#   - uses kernel_h_sigma instead of kernel_h_mu
#   - takes Z_list (squared residuals) instead of Y_list
#
# Definition:
#
#   R_r(t) = sum_{i=1}^{n} w_i sum_{j=1}^{m_i}
#              K_{h_sigma}(T_ij - t) * {(T_ij - t) / h_sigma}^r * Z_ij
#
# where Z_ij = {Y_ij - mu_hat(T_ij)}^2 are the squared residuals.
#
# Dependencies:
#   - compute_weights()   defined in weights.R
#   - kernel_h_sigma()    defined in kernel_h_sigma.R
# =============================================================================


#' Compute R_r(t) for varsigma^2 estimation at one or more evaluation points
#'
#' @param t_vec    Numeric vector. Points at which R_r is evaluated.
#'                 Can be a scalar for single-point evaluation (backward compatible).
#' @param T_list   List of length n. T_list((i)) is the vector of observation
#'                 times for subject i.
#' @param Z_list   List of length n. Z_list((i)) is the vector of squared
#'                 residuals Z_ij = {Y_ij - mu_hat(T_ij)}^2 for subject i.
#' @param h_sigma  Numeric scalar. Bandwidth h_sigma > 0.
#' @param r        Integer. Order of R_r. Must be 0, 1, or 2.
#' @param kernel   Character string. Kernel type, passed to kernel_h_sigma().
#' @param scheme   Character string. Weighting scheme passed to compute_weights().
#'                 One of "OBS", "SUBJ", "OPTIMAL".
#'
#' @return Numeric vector of length(t_vec). Returns a scalar if t_vec is scalar.

#' @export
compute_Rr_sigma <- function(t_vec, T_list, Z_list, h_sigma, r, kernel, scheme) {
  
  # --- input checks -----------------------------------------------------------
  if (!is.numeric(t_vec) || length(t_vec) < 1)
    stop("t_vec must be a non-empty numeric vector.")
  if (!is.list(T_list))
    stop("T_list must be a list.")
  if (!is.list(Z_list))
    stop("Z_list must be a list.")
  if (length(T_list) != length(Z_list))
    stop("T_list and Z_list must have the same length n.")
  if (!is.numeric(h_sigma) || length(h_sigma) != 1 || h_sigma <= 0)
    stop("h_sigma must be a single strictly positive number.")
  if (!r %in% c(0, 1, 2))
    stop("r must be 0, 1, or 2.")
  
  # --- flatten all observations and residuals ---------------------------------
  m_vec  <- sapply(T_list, length)
  T_all  <- unlist(T_list, use.names = FALSE)       # length N = sum(m_i)
  Z_all  <- unlist(Z_list, use.names = FALSE)       # length N
  
  # expand subject weights to observation level, then multiply by Z_ij
  w_subj <- compute_weights(m_vec = m_vec, scheme = scheme, h = h_sigma)
  w_all  <- rep(w_subj, times = m_vec)
  wZ_all <- w_all * Z_all
  
  # --- fully vectorized computation -------------------------------------------
  # We want D[k, j] = T_all[j] - t_vec[k], matching the theoretical formula.
  D   <- outer(t_vec, T_all, function(t, x) x - t)
  U   <- D / h_sigma
  K   <- kernel_h_sigma(D, h_sigma = h_sigma, kernel = kernel)
  wKZ <- sweep(K, 2, wZ_all, "*")
  
  # R_r[k] = sum_j wKZ[k,j] * U[k,j]^r
  Rr <- if (r == 0) rowSums(wKZ) else rowSums(wKZ * U^r)
  
  # backward compatible: return scalar if input was scalar
  if (length(t_vec) == 1) return(Rr[1])
  return(Rr)
}
