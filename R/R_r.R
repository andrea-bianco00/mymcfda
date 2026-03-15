# =============================================================================
# R_r.R
# -----------------------------------------------------------------------------
# Computes the quantity R_r(t) used in the ridged local linear estimator
# of the mean function mu(t), following Lin & Wang (2022), Section 2.1.
#
# Definition:
#
#   R_r(t) = sum_{i=1}^n w_i sum_{j=1}^{m_i}
#              K_{h_mu}(T_ij - t) * ((T_ij - t)/h_mu)^r * Y_ij
#
# for r = 0, 1.
#
# This implementation is vectorized over t_vec:
#   - input t_vec can contain one or many evaluation points
#   - output is a numeric vector of length(t_vec)
#
# Dependencies:
#   - compute_weights()  defined in weights.R
#   - kernel_h_mu()      defined in kernel_h_mu.R
# =============================================================================


#' Compute R_r(t) for one or more evaluation points
#'
#' @param t_vec   Numeric vector. Evaluation points t.
#' @param T_list  List of length n. T_list[[i]] is the vector of observation
#'                times for subject i.
#' @param Y_list  List of length n. Y_list[[i]] is the vector of observed values
#'                for subject i.
#' @param h_mu    Numeric scalar. Bandwidth h_mu > 0.
#' @param r       Integer. Order of the moment, typically 0 or 1.
#' @param kernel  Character string. Kernel type passed to kernel_h_mu().
#' @param scheme  Character string. Weighting scheme passed to compute_weights().
#'                One of "OBS", "SUBJ", "OPTIMAL".
#'
#' @return Numeric vector of length(t_vec), where the k-th entry is R_r(t_vec[k]).

compute_Rr <- function(t_vec, T_list, Y_list, h_mu, r, kernel, scheme) {
  
  # --- input checks -----------------------------------------------------------
  if (!is.numeric(t_vec) || length(t_vec) < 1)
    stop("t_vec must be a non-empty numeric vector.")
  if (!is.list(T_list) || !is.list(Y_list))
    stop("T_list and Y_list must be lists.")
  if (length(T_list) != length(Y_list))
    stop("T_list and Y_list must have the same length n.")
  if (!is.numeric(h_mu) || length(h_mu) != 1 || h_mu <= 0)
    stop("h_mu must be a single strictly positive number.")
  if (!is.numeric(r) || length(r) != 1 || r < 0)
    stop("r must be a non-negative integer.")
  
  # --- basic quantities -------------------------------------------------------
  m_vec <- sapply(T_list, length)
  w_vec <- compute_weights(m_vec = m_vec, scheme = scheme, h = h_mu)
  
  # flatten subject-specific times/values and attach subject weights
  T_all <- unlist(T_list, use.names = FALSE)
  Y_all <- unlist(Y_list, use.names = FALSE)
  W_all <- rep(w_vec, times = m_vec)
  
  # --- build matrix of shifted/scaled distances -------------------------------
  # We want D[k, j] = T_all[j] - t_vec[k], exactly as in the theoretical formula.
  D <- outer(t_vec, T_all, function(t, x) x - t)
  U <- D / h_mu
  
  # kernel matrix: K_h(T_ij - t)
  K_mat <- kernel_h_mu(D, h_mu = h_mu, kernel = kernel)
  
  # --- compute R_r(t) for each t in t_vec ------------------------------------
  R_r <- rowSums(K_mat * (U^r) *
                   matrix(W_all * Y_all,
                          nrow = length(t_vec),
                          ncol = length(T_all),
                          byrow = TRUE))
  
  # backward compatible: return scalar if input was scalar
  if (length(t_vec) == 1) return(R_r[1])
  return(R_r)
}
