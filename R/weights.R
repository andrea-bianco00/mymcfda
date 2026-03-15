# =============================================================================
# weights.R
# -----------------------------------------------------------------------------
# Weighting schemes for mean function estimation in functional data analysis.
#
# Reference: Zhang & Wang (2018) - "Optimal weighting schemes for longitudinal
#            and functional data", Statistics and Probability Letters.
#
# Three schemes are implemented:
#
#   OBS  (equal-weight-per-observation):
#        w_i = 1 / sum(m_i)
#        All observations receive the same weight regardless of subject.
#        Asymptotically optimal for non-dense data (m_i = o(n^{1/4})).
#        See Corollary 3.1(i) in Zhang & Wang (2018).
#
#   SUBJ (equal-weight-per-subject):
#        w_i = 1 / (n * m_i)
#        All subjects receive the same total weight 1/n.
#        Asymptotically optimal for ultra-dense data (m_i >> n^{1/4}).
#        See Corollary 3.1(ii) in Zhang & Wang (2018).
#
#   OPTIMAL:
#        w_i = {h^{-1} + (m_i - 1)}^{-1} / sum_j m_j {h^{-1} + (m_j-1)}^{-1}
#        Minimizes the L2 rate of convergence for any design.
#        Requires bandwidth h as input since weights depend on h.
#        See Theorem 3.1 in Zhang & Wang (2018).
#
# Normalisation constraint: sum_i m_i * w_i = 1  (required by local linear
# smoothing theory, see equation (1) in Zhang & Wang (2018)).
# =============================================================================


#' Compute weighting scheme for mean function estimation
#'
#' @param m_vec  Integer vector of length n. Number of observations per subject.
#' @param scheme Character string. One of "OBS", "SUBJ", "OPTIMAL".
#' @param h      Numeric scalar. Bandwidth h_mu. Required only for "OPTIMAL".
#'               Ignored for "OBS" and "SUBJ". Default is NULL.
#'
#' @return Numeric vector of length n with weights w_1, ..., w_n satisfying
#'         the normalisation constraint sum_i m_i * w_i = 1.
#'
#' @examples
#' m_vec <- c(3, 4, 5, 2, 4)
#' compute_weights(m_vec, scheme = "OBS")
#' compute_weights(m_vec, scheme = "SUBJ")
#' compute_weights(m_vec, scheme = "OPTIMAL", h = 0.3)

compute_weights <- function(m_vec, scheme, h = NULL) {
  
  # --- input checks -----------------------------------------------------------
  if (!is.numeric(m_vec) || any(m_vec < 1)) {
    stop("m_vec must be a numeric vector with all entries >= 1.")
  }
  if (!scheme %in% c("OBS", "SUBJ", "OPTIMAL")) {
    stop("scheme must be one of: 'OBS', 'SUBJ', 'OPTIMAL'.")
  }
  if (scheme == "OPTIMAL" && is.null(h)) {
    stop("h must be provided when scheme = 'OPTIMAL'.")
  }
  if (!is.null(h) && (h <= 0)) {
    stop("h must be strictly positive.")
  }
  
  # --- compute weights --------------------------------------------------------
  n <- length(m_vec)
  
  if (scheme == "OBS") {
    # -------------------------------------------------------------------------
    # OBS: equal weight per observation
    # w_i = 1 / sum_j m_j
    # Normalisation: sum_i m_i * w_i = sum_i m_i / sum_j m_j = 1  [OK]
    # -------------------------------------------------------------------------
    w <- rep(1 / sum(m_vec), n)
    
  } else if (scheme == "SUBJ") {
    # -------------------------------------------------------------------------
    # SUBJ: equal weight per subject
    # w_i = 1 / (n * m_i)
    # Normalisation: sum_i m_i * (1 / (n * m_i)) = n / n = 1  [OK]
    # -------------------------------------------------------------------------
    w <- 1 / (n * m_vec)
    
  } else if (scheme == "OPTIMAL") {
    # -------------------------------------------------------------------------
    # OPTIMAL: minimises L2 rate of convergence (Theorem 3.1, Zhang & Wang 2018)
    # w_i = {h^{-1} + (m_i - 1)}^{-1} / sum_j m_j {h^{-1} + (m_j - 1)}^{-1}
    # -------------------------------------------------------------------------
    raw   <- 1 / (1/h + (m_vec - 1))          # {h^{-1} + (m_i - 1)}^{-1}
    denom <- sum(m_vec * raw)                  # sum_j m_j * raw_j
    w     <- raw / denom
  }
  
  # --- verify normalisation constraint ----------------------------------------
  # sum_i m_i * w_i should equal 1 (up to numerical precision)
  constraint <- sum(m_vec * w)
  if (abs(constraint - 1) > 1e-10) {
    warning(sprintf(
      "Normalisation constraint violated: sum(m_i * w_i) = %.6f (expected 1).",
      constraint
    ))
  }
  
  return(w)
}
