# =============================================================================
# compute_residuals.R
# -----------------------------------------------------------------------------
# Computes the squared residuals Z_ij = {Y_ij - mu_hat(T_ij)}^2 for each
# subject i and observation j.
#
# These residuals are the input for the local linear smoother of varsigma^2(t)
# in Block 2, as described in Section 2.2 of Lin & Wang (2022).
#
# Note: mu_hat is evaluated at the observed time points T_ij, NOT on a grid.
#
# OPTIMIZED IMPLEMENTATION:
#   Instead of calling compute_mu_hat_t() once per observation (n * m_bar
#   calls total), all observed time points are flattened into a single vector
#   T_all of length N = sum(m_i), and compute_mu_hat_t() is called once:
#
#     T_all    <- unlist(T_list)             # all N observation times
#     mu_all   <- compute_mu_hat_t(T_all, .) # mu_hat at all N points at once
#     Z_all    <- (Y_all - mu_all)^2         # squared residuals, length N
#     Z_list   <- relist(Z_all, T_list)      # restore list structure
#
#   This reduces n * m_bar kernel evaluations to a single matrix operation.
#
# Dependencies:
#   - compute_mu_hat_t()  defined in mu_hat_t.R
#   - compute_Sr()        defined in S_r.R
#   - compute_Rr()        defined in R_r.R
#   - compute_weights()   defined in weights.R
#   - kernel_h_mu()       defined in kernel_h_mu.R
#   - kernel_fun()        defined in kernel.R
# =============================================================================


#' Compute squared residuals Z_ij = {Y_ij - mu_hat(T_ij)}^2
#'
#' @param T_list  List of length n. T_list[[i]] is the vector of observation
#'                times for subject i.
#' @param Y_list  List of length n. Y_list[[i]] is the vector of observed
#'                values Y_ij for subject i.
#' @param h_mu    Numeric scalar. Bandwidth used to estimate mu_hat.
#' @param kernel  Character string. Kernel type, passed to compute_mu_hat_t().
#' @param scheme  Character string. Weighting scheme passed to compute_mu_hat_t().
#'                One of "OBS", "SUBJ", "OPTIMAL".
#'
#' @return A list of length n. Z_list[[i]] is the vector of squared residuals
#'         Z_ij = {Y_ij - mu_hat(T_ij)}^2 for subject i.

compute_residuals <- function(T_list, Y_list, h_mu, kernel, scheme) {
  
  # --- input checks -----------------------------------------------------------
  if (!is.list(T_list) || !is.list(Y_list))
    stop("T_list and Y_list must be lists.")
  if (length(T_list) != length(Y_list))
    stop("T_list and Y_list must have the same length n.")
  if (!is.numeric(h_mu) || length(h_mu) != 1 || h_mu <= 0)
    stop("h_mu must be a strictly positive scalar.")
  
  # --- flatten all observed points and responses ------------------------------
  m_vec <- sapply(T_list, length)
  T_all <- unlist(T_list)          # all N = sum(m_i) observation times
  Y_all <- unlist(Y_list)          # corresponding responses
  
  # --- evaluate mu_hat at all observed points in one vectorized call ----------
  mu_all <- compute_mu_hat_t(t_vec  = T_all,
                             T_list = T_list,
                             Y_list = Y_list,
                             h_mu   = h_mu,
                             kernel = kernel,
                             scheme = scheme)
  
  # --- squared residuals Z_ij = {Y_ij - mu_hat(T_ij)}^2 ----------------------
  Z_all <- (Y_all - mu_all)^2
  
  # --- restore list structure: Z_list[[i]] has length m_i ---------------------
  Z_list <- split(Z_all, rep(seq_along(m_vec), times = m_vec))
  names(Z_list) <- NULL
  
  return(Z_list)
}
