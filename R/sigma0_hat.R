# =============================================================================
# sigma0_hat.R
# -----------------------------------------------------------------------------
# Computes the noise variance estimator sigma^2_0 of Lin & Wang (2022),
# Section 3, using the ridged version of equation (8):
#
#   sigma^2_0_hat = (A_hat_0 - A_hat_2) / (B_hat + Delta)
#
# where Delta = (n * m_bar)^{-2} * h0 is a ridging term that prevents
# division by zero when B_hat is very small.
#
# Dependencies:
#   - compute_A0()  defined in A0_hat.R
#   - compute_A2()  defined in A2_hat.R
#   - compute_B()   defined in B_hat.R
# =============================================================================

compute_sigma0 <- function(T_list, Y_list, h0) {
  
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
  
  # --- basic quantities -------------------------------------------------------
  m_vec <- sapply(T_list, length)
  n     <- length(m_vec)
  m_bar <- mean(m_vec)
  
  # --- ridging term -----------------------------------------------------------
  Delta <- (n * m_bar)^(-2) * h0
  
  # --- compute A_hat_0, A_hat_2, B_hat ---------------------------------------
  A0 <- compute_A0(T_list, Y_list, h0)
  A2 <- compute_A2(T_list, Y_list, h0)
  B  <- compute_B(T_list, h0)
  
  # --- compute sigma^2_0_hat --------------------------------------------------
  sigma0_hat <- (A0 - A2) / (B + Delta)
  
  cat(sprintf("A_hat_0    = %.6f\n", A0))
  cat(sprintf("A_hat_2    = %.6f\n", A2))
  cat(sprintf("A0 - A2    = %.6f\n", A0 - A2))
  cat(sprintf("B_hat      = %.6f\n", B))
  cat(sprintf("Delta      = %.2e\n", Delta))
  cat(sprintf("sigma^2_0  = %.6f\n", sigma0_hat))
  
  return(list(
    sigma0_hat = sigma0_hat,
    A0         = A0,
    A2         = A2,
    B          = B,
    Delta      = Delta
  ))
}
