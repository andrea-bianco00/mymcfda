# =============================================================================
# draw_measurement_times.R
# -----------------------------------------------------------------------------
# Implements Step 3 of the data-generating mechanism of Lin & Wang (2022),
# Section 1.7:
#
#   "Conditional on O_i, draw T_i1, ..., T_im_i i.i.d. from the density
#    f_{T|O}(t | O_i) = f_0(t - O_i + delta/2),
#    where f_0 > 0 on [0, delta] and f_0' is Lipschitz continuous on [0, delta]."
#
# SHIFT STRUCTURE (Assumption 3, Section 1.6):
#   f_0 is a template density on the reference interval [0, delta].
#   Given O_i, the measurement times are obtained by shifting f_0 to the
#   subject's window [O_i - delta/2, O_i + delta/2]:
#
#     S_j  ~ f_0 on [0, delta]             (relative position within window)
#     T_ij  = S_j + O_i - delta/2          (absolute position on [0, 1])
#
#   so that T_ij in [O_i - delta/2, O_i + delta/2] subset [0, 1].
#
# This function does NOT contain any mathematical or sampling logic.
# It delegates the sampling of S_j entirely to sample_fzero.R, which in
# turn calls the appropriate density-specific file. The shift is applied here.
#
# INPUT:
#   O_vec  : numeric vector of length n with reference times O_1,...,O_n.
#             Typically produced by draw_reference_times().
#   m_vec  : integer vector of length n with the number of measurements
#             m_i >= 2 for each subject i. m_i >= 2 is required to enable
#             estimation of sigma_0^2 via within-subject pairs (Section 3.2.2).
#   delta  : snippet length in (0, 1).
#   type   : family of f_0 (passed to sample_fzero).
#   params : named list of parameters for the chosen family (passed to
#            sample_fzero). See sample_fzero.R for full documentation.
#
# OUTPUT:
#   A list of n numeric vectors. The i-th element T_list[[i]] is a vector
#   of length m_i containing the measurement times T_i1, ..., T_im_i for
#   subject i, all in [O_i - delta/2, O_i + delta/2] subset [0, 1].
#   These times ARE observed by the analyst (unlike O_i).
#
# REQUIRED SOURCES:
#   source("sample_fzero.R")   which sources all sample_fzero_*.R files
#
# Dependencies: base R only (via the sourced files).
# =============================================================================



#' Draw measurement times T_i1,...,T_im_i for each subject (Step 3 of Lin & Wang 2022)
#'
#' @param O_vec  Numeric vector of length n. Reference times O_1,...,O_n,
#'               all in (delta/2, 1 - delta/2). Produced by
#'               draw_reference_times().
#' @param m_vec  Integer vector of length n. Number of measurements m_i >= 2
#'               for each subject i. m_i >= 2 is required for estimation of
#'               sigma_0^2 via within-subject pairs (Section 3.2.2 of the paper).
#'               Can be a scalar (same m for all subjects) or a vector.
#' @param delta  Numeric in (0, 1). Snippet length. Determines the window
#'               (O_i - delta/2, O_i + delta/2) for each subject.
#' @param type   Character. Family of f_0. One of:
#'                 "uniform"   (default, as in Lin & Wang 2022)
#'                 "truncnorm"
#'                 "truncbeta"
#'                 "truncexp"
#' @param params Named list of parameters for the chosen f_0 family.
#'               Passed directly to sample_fzero(). See sample_fzero.R for
#'               the full parameter documentation of each family.
#'
#' @return A list of n numeric vectors. T_list((i)) contains the m_i
#'         measurement times for subject i, all in
#'         (O_i - delta/2, O_i + delta/2) subset (0, 1).
#'         These times are OBSERVED by the analyst.
#'
#' @examples
#' # Step 2: draw reference times
#' ref  <- draw_reference_times(n = 5, delta = 0.4)
#'
#' # Step 3a: all subjects have the same number of measurements
#' T_list <- draw_measurement_times(O_vec = ref$O_vec, m_vec = 4, delta = 0.4)
#' lengths(T_list)   # all 4
#'
#' # Step 3b: heterogeneous measurement counts
#' T_list2 <- draw_measurement_times(O_vec  = ref$O_vec,
#'                                   m_vec  = c(2, 3, 4, 5, 6),
#'                                   delta  = 0.4)
#' lengths(T_list2)  # 2 3 4 5 6
#'
#' # Check that all T_ij lie inside their subject's window
#' for (i in seq_along(T_list2)) {
#'   A_i <- ref$O_vec[i] - 0.4 / 2
#'   B_i <- ref$O_vec[i] + 0.4 / 2
#'   stopifnot(all(T_list2[[i]] >= A_i - 1e-10 & T_list2[[i]] <= B_i + 1e-10))
#' }
#'
#' # Truncated Normal f_0
#' T_list3 <- draw_measurement_times(O_vec  = ref$O_vec,
#'                                   m_vec  = 4,
#'                                   delta  = 0.4,
#'                                   type   = "truncnorm",
#'                                   params = list(mu_0 = 0.2, sd_0 = 0.08))

#' @export
draw_measurement_times <- function(O_vec,
                                   m_vec,
                                   delta,
                                   type   = "uniform",
                                   params = list()) {

  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.numeric(O_vec) || length(O_vec) < 1)
    stop("O_vec must be a non-empty numeric vector.")
  if (!is.numeric(delta) || length(delta) != 1 || delta <= 0 || delta >= 1)
    stop("delta must be a number in (0, 1).")

  n <- length(O_vec)

  # expand scalar m_vec to a vector of length n
  if (length(m_vec) == 1) {
    m_vec <- rep(as.integer(m_vec), n)
  } else {
    m_vec <- as.integer(m_vec)
  }

  if (length(m_vec) != n)
    stop("m_vec must be a scalar or a vector of the same length as O_vec.")
  if (any(m_vec < 2))
    stop("All m_i must be >= 2 (required for within-subject pair estimation ",
         "of sigma_0^2 in Section 3.2.2 of Lin & Wang 2022).")

  # check O_vec values are in [delta/2, 1 - delta/2]
  lo <- delta / 2
  hi <- 1 - delta / 2
  if (any(O_vec < lo - 1e-10) || any(O_vec > hi + 1e-10))
    stop("All O_i must lie in [delta/2, 1 - delta/2]. ",
         "Use draw_reference_times() to generate valid O_vec.")

  if (!type %in% c("uniform", "truncnorm", "truncbeta", "truncexp"))
    stop("type must be one of: 'uniform', 'truncnorm', 'truncbeta', 'truncexp'.")
  if (!is.list(params))
    stop("params must be a named list.")

  # ---------------------------------------------------------------------------
  # STEP 3 (Lin & Wang 2022, Section 1.7):
  # For each subject i:
  #   1. Draw S_j ~ f_0 on [0, delta]    (relative positions)
  #   2. Shift: T_ij = S_j + O_i - delta/2  (absolute positions on [0, 1])
  # ---------------------------------------------------------------------------
  T_list <- vector("list", n)

  for (i in seq_len(n)) {

    # draw m_i relative positions from f_0 on [0, delta]
    S_i <- sample_fzero(m      = m_vec[i],
                        delta  = delta,
                        type   = type,
                        params = params)

    # shift to subject i's window [O_i - delta/2, O_i + delta/2]
    T_list[[i]] <- S_i + O_vec[i] - delta / 2
  }

  # ---------------------------------------------------------------------------
  # SANITY CHECK: all T_ij must lie inside [0, 1]
  # ---------------------------------------------------------------------------
  all_T <- unlist(T_list)
  if (any(all_T < -1e-10) || any(all_T > 1 + 1e-10))
    stop("Internal error: some measurement times fall outside [0, 1].")

  return(T_list)
}
