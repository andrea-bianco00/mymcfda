# =============================================================================
# varsigma2_hat_t.R
# -----------------------------------------------------------------------------
# Computes the ridged local linear estimator of varsigma^2(t) at one or more
# points.
#
# varsigma^2(t) = E{Y(t) - mu(t)}^2 = sigma^2_X(t) + sigma^2_0
#
# This is the Block 2 analogue of mu_hat_t.R. The formula is identical to
# the ridged local linear estimator in equation (2) of Lin & Wang (2022),
# but applied to the squared residuals Z_ij = {Y_ij - mu_hat(T_ij)}^2
# with bandwidth h_sigma instead of h_mu.
#
# Definition:
#
#   varsigma2_hat(t) = (R_0*S_2 - R_1*S_1)
#                      / (S_0*S_2 - S_1^2 + Delta * 1_{|S_0*S_2-S_1^2|<Delta})
#
# where:
#   S_r(t) = sum_i w_i sum_j K_{h_sigma}(T_ij-t) * {(T_ij-t)/h_sigma}^r
#   R_r(t) = sum_i w_i sum_j K_{h_sigma}(T_ij-t) * {(T_ij-t)/h_sigma}^r * Z_ij
#   Delta  = (n * m_bar)^{-2}
#
# OPTIMIZED IMPLEMENTATION:
#   Accepts a vector t_vec of evaluation points. Each of S_r_sigma and
#   R_r_sigma is computed in a single vectorized call (returning a vector of
#   length n_t), and the final formula is applied element-wise.
#   Backward compatible: scalar t returns scalar varsigma2_hat.
#
# Dependencies:
#   - compute_Sr_sigma()  defined in S_r_sigma.R
#   - compute_Rr_sigma()  defined in R_r_sigma.R
#   - compute_weights()   defined in weights.R
#   - kernel_h_sigma()    defined in kernel_h_sigma.R
#   - kernel_fun()        defined in kernel.R
# =============================================================================


#' Ridged local linear estimator of varsigma^2(t) at one or more points
#'
#' @param t_vec    Numeric vector. Points at which varsigma^2(t) is estimated.
#'                 Can be a scalar for single-point evaluation (backward compatible).
#' @param T_list   List of length n. T_list((i)) is the vector of observation
#'                 times for subject i.
#' @param Z_list   List of length n. Z_list((i)) is the vector of squared
#'                 residuals Z_ij = {Y_ij - mu_hat(T_ij)}^2 for subject i.
#' @param h_sigma  Numeric scalar. Bandwidth h_sigma > 0.
#' @param kernel   Character string. Kernel type, passed to kernel_h_sigma().
#' @param scheme   Character string. Weighting scheme passed to compute_weights().
#'                 One of "OBS", "SUBJ", "OPTIMAL".
#'
#' @return Numeric vector of length(t_vec). Returns a scalar if t_vec is scalar.

#' @export
compute_varsigma2_hat_t <- function(t_vec, T_list, Z_list, h_sigma, kernel, scheme) {
  
  # --- input checks -----------------------------------------------------------
  if (!is.numeric(t_vec) || length(t_vec) < 1)
    stop("t_vec must be a non-empty numeric vector.")
  if (!is.list(T_list) || !is.list(Z_list))
    stop("T_list and Z_list must be lists.")
  if (length(T_list) != length(Z_list))
    stop("T_list and Z_list must have the same length n.")
  
  # --- Delta = (n * m_bar)^{-2} -----------------------------------------------
  m_vec <- sapply(T_list, length)
  Delta <- (length(m_vec) * mean(m_vec))^(-2)
  
  # --- S_0, S_1, S_2 and R_0, R_1 — each one vectorized call -----------------
  S0 <- compute_Sr_sigma(t_vec = t_vec, T_list = T_list,
                         h_sigma = h_sigma, r = 0, kernel = kernel, scheme = scheme)
  S1 <- compute_Sr_sigma(t_vec = t_vec, T_list = T_list,
                         h_sigma = h_sigma, r = 1, kernel = kernel, scheme = scheme)
  S2 <- compute_Sr_sigma(t_vec = t_vec, T_list = T_list,
                         h_sigma = h_sigma, r = 2, kernel = kernel, scheme = scheme)
  
  R0 <- compute_Rr_sigma(t_vec = t_vec, T_list = T_list, Z_list = Z_list,
                         h_sigma = h_sigma, r = 0, kernel = kernel, scheme = scheme)
  R1 <- compute_Rr_sigma(t_vec = t_vec, T_list = T_list, Z_list = Z_list,
                         h_sigma = h_sigma, r = 1, kernel = kernel, scheme = scheme)
  
  # --- ridged local linear estimator ------------------------------------------
  numerator   <- R0 * S2 - R1 * S1
  denom_raw   <- S0 * S2 - S1^2
  denominator <- denom_raw + Delta * as.numeric(abs(denom_raw) < Delta)
  
  varsigma2_hat <- numerator / denominator
  
  # backward compatible: return scalar if input was scalar
  if (length(t_vec) == 1) return(varsigma2_hat[1])
  return(varsigma2_hat)
}
