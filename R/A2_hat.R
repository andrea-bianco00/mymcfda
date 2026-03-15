# ============================================================
# compute_A2
# ------------------------------------------------------------
# Compute the quantity A2_hat used in the estimation of sigma0^2.
#
# Formula (Lin & Wang):
#
# A2_hat =
# (1/n) * sum_i [ 1 / (m_i (m_i - 1)) *
#                 sum_{j != l} Y_ij * Y_il * 1{|T_ij - T_il| < h0} ]
#
# Inputs
#   T_list : list of observation times per subject
#   Y_list : list of responses per subject
#   h0     : neighbourhood radius
#
# Output
#   scalar A2_hat
# ============================================================

compute_A2 <- function(T_list, Y_list, h0) {
  
  n <- length(T_list)
  
  total <- 0
  
  for (i in seq_len(n)) {
    
    Ti <- T_list[[i]]
    Yi <- Y_list[[i]]
    
    mi <- length(Ti)
    
    if (mi < 2) next
    
    sum_pairs <- 0
    
    for (j in seq_len(mi)) {
      for (l in seq_len(mi)) {
        
        if (j == l) next
        
        if (abs(Ti[j] - Ti[l]) < h0) {
          sum_pairs <- sum_pairs + Yi[j] * Yi[l]
        }
        
      }
    }
    
    total <- total + (1 / (mi * (mi - 1))) * sum_pairs
  }
  
  A2_hat <- total / n
  
  return(A2_hat)
}
