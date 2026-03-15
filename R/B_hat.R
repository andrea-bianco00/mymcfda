# =============================================================================
# B_hat.R
# -----------------------------------------------------------------------------
# Computes the quantity B_hat used in the noise variance estimator
# sigma^2_0 of Lin & Wang (2022), Section 3.
#
# Definition (equation 6):
#
#   B_hat = (1/n) * sum_{i=1}^{n} (1 / (m_i*(m_i-1))) *
#             sum_{j != l} 1_{|T_ij - T_il| < h0}
#
# Intuition: B_hat estimates
#
#   B = E[1_{|T1 - T2| < h0}]
#
# which is the probability that two observation times from the same subject
# are within h0 of each other. It acts as the normalizing constant in the
# estimator sigma^2_0 = (A_hat_0 - A_hat_1) / B_hat.
#
# B_hat is always in [0, 1] since it is an average of indicators.
#
# Dependencies: none (only base R)
# =============================================================================


#' Compute B_hat
#'
#' @param T_list  List of length n. T_list((i)) is the vector of observation
#'                times for subject i.
#' @param h0      Numeric scalar. Bandwidth h0 > 0, computed via compute_h0().
#'
#' @return Numeric scalar. The value of B_hat, always in (0, 1).
#'
#' @examples
#' T_list <- list(c(0.1, 0.3, 0.4), c(0.5, 0.6, 0.7, 0.8), c(0.2, 0.4))
#' compute_B(T_list, h0 = 0.15)

#' @export
compute_B <- function(T_list, h0) {
  
  # --- input checks -----------------------------------------------------------
  if (!is.list(T_list)) {
    stop("T_list must be a list.")
  }
  if (!is.numeric(h0) || length(h0) != 1 || h0 <= 0) {
    stop("h0 must be a strictly positive scalar.")
  }
  
  n <- length(T_list)
  
  # --- compute B_hat ----------------------------------------------------------
  # B_hat = (1/n) * sum_i (1/(m_i*(m_i-1))) * sum_{j!=l} 1_{|T_ij-T_il|<h0}
  B <- 0
  
  for (i in seq_len(n)) {
    T_i <- T_list[[i]]
    m_i <- length(T_i)
    
    if (m_i < 2) next
    
    # count ordered pairs (j, l) with j != l and |T_ij - T_il| < h0
    count <- 0
    for (j in seq_len(m_i)) {
      for (l in seq_len(m_i)) {
        if (j != l && abs(T_i[j] - T_i[l]) < h0) {
          count <- count + 1
        }
      }
    }
    
    B <- B + count / (m_i * (m_i - 1))
  }
  
  return(B / n)
}
