# =============================================================================
# h0_rule.R
# -----------------------------------------------------------------------------
# Computes the bandwidth h0 for the noise variance estimator sigma^2_0
# using the empirical rule described in Section 3 and Appendix A.1 of
# Lin & Wang (2022):
#
#   h0 = 0.29 * delta_hat * ||varsigma2_hat||_2 * (n * m_bar^2)^{-1/5}
#
# where:
#   - delta_hat         = max_i max_{j != l} |T_ij - T_il|
#                         (estimated snippet width)
#   - ||varsigma2_hat||_2 = integral of varsigma2_hat(t) dt
#                           (overall variability of the data)
#   - (n * m_bar^2)^{-1/5} = theoretically optimal rate from Theorem 1
#
# After computing h0, a safety check is performed: if the neighborhood
#
#   N(h0) = {(T_ij, T_il) : |T_ij - T_il| < h0, i=1,...,n, j != l}
#
# is empty or contains fewer than floor(0.1 * sum_i m_i*(m_i-1)) pairs,
# h0 is increased until the condition is satisfied, as recommended in
# Section 3 of Lin & Wang (2022).
#
# Dependencies: none (only base R)
# =============================================================================


#' Compute h0 via the empirical rule of Lin & Wang (2022)
#'
#' @param T_list        List of length n. T_list((i)) is the vector of
#'                      observation times for subject i.
#' @param res_varsigma  Output of compute_varsigma2_hat(). A list with
#'                      components $t_grid and $varsigma2_hat.
#'
#' @return A list with components:
#'   - h0:          Numeric scalar. The selected bandwidth.
#'   - delta_hat:   Numeric scalar. Estimated snippet width.
#'   - norm_varsig: Numeric scalar. ||varsigma2_hat||_2.
#'   - n_pairs:     Integer. Number of pairs in N(h0).
#'   - adjusted:    Logical. TRUE if h0 was increased by the safety check.
#'
#' @examples
#' T_list <- list(c(0.1, 0.3, 0.4), c(0.5, 0.6, 0.7, 0.8), c(0.2, 0.4))
#' # assume res_varsigma is the output of compute_varsigma2_hat()
#' # h0_res <- compute_h0(T_list, res_varsigma)

#' @export
compute_h0 <- function(T_list, res_varsigma) {

  # --- input checks -----------------------------------------------------------
  if (!is.list(T_list)) {
    stop("T_list must be a list.")
  }
  if (!is.list(res_varsigma) ||
      !all(c("t_grid", "varsigma2_hat") %in% names(res_varsigma))) {
    stop("res_varsigma must be a list with components $t_grid and $varsigma2_hat.")
  }

  # --- basic quantities -------------------------------------------------------
  m_vec <- sapply(T_list, length)
  n     <- length(m_vec)
  m_bar <- mean(m_vec)

  # --- Step 1: compute delta_hat ----------------------------------------------
  # delta_hat = max_i max_{j != l} |T_ij - T_il|
  # = max_i (max(T_i) - min(T_i))
  delta_hat <- max(sapply(T_list, function(t_i) max(t_i) - min(t_i)))

  # --- Step 2: compute ||varsigma2_hat||_2 ------------------------------------
  # For the paper: ||sigma_hat||^2_2 = integral sigma_hat^2(t) dt  (p. 6, Lin & Wang 2022)
  # therefore      ||sigma_hat||_2   = sqrt(integral sigma_hat^2(t) dt)
  # Note: varsigma2_hat estimates sigma^2(t), so sigma_hat(t) = sqrt(varsigma2_hat(t))
  # and integral sigma_hat^2(t) dt = integral varsigma2_hat(t) dt
  # The correct norm is therefore: sqrt(integral varsigma2_hat(t) dt)
  t_grid        <- res_varsigma$t_grid
  varsigma2_hat <- res_varsigma$varsigma2_hat
  norm_varsig   <- sqrt(sum((varsigma2_hat[-1] + varsigma2_hat[-length(varsigma2_hat)]) / 2 *
                              diff(t_grid)))

  # --- Step 3: compute h0 with empirical rule ---------------------------------
  h0 <- 0.29 * delta_hat * norm_varsig * (n * m_bar^2)^(-1/5)

  cat(sprintf("delta_hat       = %.4f\n", delta_hat))
  cat(sprintf("||varsigma2||_2 = %.4f\n", norm_varsig))
  cat(sprintf("(n*m_bar^2)^(-1/5) = %.4f\n", (n * m_bar^2)^(-1/5)))
  cat(sprintf("h0 (empirical rule) = %.4f\n", h0))

  # --- Step 4: safety check ---------------------------------------------------
  # Minimum number of pairs required in N(h0)
  total_pairs  <- sum(m_vec * (m_vec - 1))
  min_pairs    <- floor(0.1 * total_pairs)

  # count pairs in N(h0)
  count_pairs <- function(h) {
    count <- 0
    for (i in seq_len(n)) {
      T_i <- T_list[[i]]
      m_i <- length(T_i)
      if (m_i < 2) next
      for (j in seq_len(m_i)) {
        for (l in seq_len(m_i)) {
          if (j != l && abs(T_i[j] - T_i[l]) < h) {
            count <- count + 1
          }
        }
      }
    }
    return(count)
  }

  n_pairs  <- count_pairs(h0)
  adjusted <- FALSE

  if (n_pairs < min_pairs) {
    cat(sprintf("Safety check: N(h0) has %d pairs < min %d - increasing h0\n",
                n_pairs, min_pairs))

    # increase h0 by 10% steps until condition is met
    h0_try <- h0
    while (n_pairs < min_pairs) {
      h0_try  <- h0_try * 1.10
      n_pairs <- count_pairs(h0_try)
    }
    h0       <- h0_try
    adjusted <- TRUE
    cat(sprintf("h0 adjusted to %.4f  (N(h0) = %d pairs)\n", h0, n_pairs))
  } else {
    cat(sprintf("Safety check passed: N(h0) = %d pairs >= min %d\n",
                n_pairs, min_pairs))
  }

  return(list(
    h0          = h0,
    delta_hat   = delta_hat,
    norm_varsig = norm_varsig,
    n_pairs     = n_pairs,
    adjusted    = adjusted
  ))
}
