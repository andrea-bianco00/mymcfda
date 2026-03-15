# =============================================================================
# A0_hat.R
# -----------------------------------------------------------------------------
# Computes the quantity A_hat_0 used in the noise variance estimator
# sigma^2_0 of Lin & Wang (2022), Section 3.
#
# Definition (equation 5):
#
#   A_hat_0 = (1/n) * sum_{i=1}^{n} (1 / (m_i*(m_i-1))) *
#               sum_{j != l} Y^2_ij * 1_{|T_ij - T_il| < h0}
#
# Intuition: A_hat_0 estimates
#
#   A_0 = E[{C(T1,T1) + mu(T1)^2 + sigma^2_0} * 1_{|T1-T2| < h0}]
#
# which includes the noise variance sigma^2_0 because E[Y^2_ij] =
# E[X_i(T_ij)^2] + sigma^2_0.
#
# Dependencies: none (only base R)
# =============================================================================


#' Compute A_hat_0
#'
#' @param T_list  List of length n. T_list((i)) is the vector of observation
#'                times for subject i.
#' @param Y_list  List of length n. Y_list((i)) is the vector of observed
#'                values Y_ij for subject i.
#' @param h0      Numeric scalar. Bandwidth h0 > 0, computed via compute_h0().
#'
#' @return Numeric scalar. The value of A_hat_0.
#'
#' @examples
#' T_list <- list(c(0.1, 0.3, 0.4), c(0.5, 0.6, 0.7, 0.8), c(0.2, 0.4))
#' Y_list <- list(c(1.2, 0.8, 1.1), c(0.9, 1.3, 1.0, 1.2), c(1.1, 0.7))
#' compute_A0(T_list, Y_list, h0 = 0.1)

#' @export
compute_A0 <- function(T_list, Y_list, h0) {
  
  # --- input checks -----------------------------------------------------------
  if (!is.list(T_list) || !is.list(Y_list)) {
    stop("T_list and Y_list must be lists.")
  }
  if (length(T_list) != length(Y_list)) {
    stop("T_list and Y_list must have the same length n.")
  }
  if (!is.numeric(h0) || length(h0) != 1 || h0 <= 0) {
    stop("h0 must be a strictly positive scalar.")
  }
  
  n <- length(T_list)
  
  # --- compute A_hat_0 --------------------------------------------------------
  # A_hat_0 = (1/n) * sum_i (1/(m_i*(m_i-1))) * sum_{j!=l} Y^2_ij * 1_{|T_ij-T_il|<h0}
  A0 <- 0
  
  for (i in seq_len(n)) {
    T_i <- T_list[[i]]
    Y_i <- Y_list[[i]]
    m_i <- length(T_i)
    
    if (m_i < 2) next
    
    # sum over all ordered pairs (j, l) with j != l
    contrib <- 0
    for (j in seq_len(m_i)) {
      for (l in seq_len(m_i)) {
        if (j != l && abs(T_i[j] - T_i[l]) < h0) {
          contrib <- contrib + Y_i[j]^2
        }
      }
    }
    
    A0 <- A0 + contrib / (m_i * (m_i - 1))
  }
  
  return(A0 / n)
}
