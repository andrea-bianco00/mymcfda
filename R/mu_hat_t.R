# =============================================================================
# mu_hat_t.R
# -----------------------------------------------------------------------------
# Computes the ridged local linear estimator of the mean function mu(t)
# at one or more points.
#
# Definition (Section 2.1, Lin & Wang 2022, equation 2):
#
#   mu_hat(t) = (R_0(t) * S_2(t) - R_1(t) * S_1(t))
#               / (S_0(t) * S_2(t) - S_1(t)^2
#                  + Delta * 1_{|S_0*S_2 - S_1^2| < Delta})
#
# where:
#   S_r(t) = sum_i w_i sum_j K_{h_mu}(T_ij - t) * {(T_ij - t) / h_mu}^r
#   R_r(t) = sum_i w_i sum_j K_{h_mu}(T_ij - t) * {(T_ij - t) / h_mu}^r * Y_ij
#   Delta  = (n * m_bar)^{-2},  m_bar = mean(m_i)
#
# The ridging term stabilises the denominator when the local design matrix
# is nearly singular (few observations in the neighbourhood of t).
#
# OPTIMIZED IMPLEMENTATION:
#   Accepts a vector t_vec of evaluation points. Each of S_r and R_r is
#   computed in a single vectorized call (returning a vector of length n_t),
#   and the final formula is applied element-wise. This replaces the previous
#   pattern of calling compute_mu_hat_t() once per grid point.
#   Backward compatible: scalar t returns scalar mu_hat.
#
# Dependencies:
#   - compute_Sr()       defined in S_r.R
#   - compute_Rr()       defined in R_r.R
#   - compute_weights()  defined in weights.R
#   - kernel_h_mu()      defined in kernel_h_mu.R
#   - kernel_fun()       defined in kernel.R
# =============================================================================


#' Ridged local linear estimator of mu(t) at one or more points
#'
#' @param t_vec   Numeric vector. Points at which mu(t) is estimated.
#'                Can be a scalar for single-point evaluation (backward compatible).
#' @param T_list  List of length n. T_list((i)) is the vector of observation
#'                times for subject i.
#' @param Y_list  List of length n. Y_list((i)) is the vector of observed
#'                values Y_ij for subject i.
#' @param h_mu    Numeric scalar. Bandwidth h_mu > 0.
#' @param kernel  Character string. Kernel type, passed to kernel_h_mu().
#' @param scheme  Character string. Weighting scheme passed to compute_weights().
#'                One of "OBS", "SUBJ", "OPTIMAL".
#'
#' @return Numeric vector of length(t_vec). Returns a scalar if t_vec is scalar.

#' @export
compute_mu_hat_t <- function(t_vec, T_list, Y_list, h_mu, kernel, scheme) {
  
  # --- input checks -----------------------------------------------------------
  if (!is.numeric(t_vec) || length(t_vec) < 1)
    stop("t_vec must be a non-empty numeric vector.")
  if (!is.list(T_list) || !is.list(Y_list))
    stop("T_list and Y_list must be lists.")
  if (length(T_list) != length(Y_list))
    stop("T_list and Y_list must have the same length n.")
  
  # --- Delta = (n * m_bar)^{-2} -----------------------------------------------
  m_vec <- sapply(T_list, length)
  Delta <- (length(m_vec) * mean(m_vec))^(-2)
  
  # --- S_0, S_1, S_2 and R_0, R_1 — each one vectorized call -----------------
  S0 <- compute_Sr(t_vec = t_vec, T_list = T_list,
                   h_mu = h_mu, r = 0, kernel = kernel, scheme = scheme)
  S1 <- compute_Sr(t_vec = t_vec, T_list = T_list,
                   h_mu = h_mu, r = 1, kernel = kernel, scheme = scheme)
  S2 <- compute_Sr(t_vec = t_vec, T_list = T_list,
                   h_mu = h_mu, r = 2, kernel = kernel, scheme = scheme)
  
  R0 <- compute_Rr(t_vec = t_vec, T_list = T_list, Y_list = Y_list,
                   h_mu = h_mu, r = 0, kernel = kernel, scheme = scheme)
  R1 <- compute_Rr(t_vec = t_vec, T_list = T_list, Y_list = Y_list,
                   h_mu = h_mu, r = 1, kernel = kernel, scheme = scheme)
  
  # --- ridged local linear estimator ------------------------------------------
  numerator   <- R0 * S2 - R1 * S1
  denom_raw   <- S0 * S2 - S1^2
  denominator <- denom_raw + Delta * as.numeric(abs(denom_raw) < Delta)
  
  mu_hat <- numerator / denominator
  
  # backward compatible: return scalar if input was scalar
  if (length(t_vec) == 1) return(mu_hat[1])
  return(mu_hat)
}
