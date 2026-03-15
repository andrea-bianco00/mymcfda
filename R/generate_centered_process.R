# =============================================================================
# generate_centered_process.R
# -----------------------------------------------------------------------------
# Unified dispatcher for generating the centred process values
#
#   Z_i(T_i1), ..., Z_i(T_im_i)
#
# through one of two roads:
#
#   road = "kl"
#     -> calls generate_KL_process()
#
#   road = "decomposition"
#     -> calls generate_decomposition_process()
#
# PURPOSE:
#   This function provides a single entry point to generate the centred
#   process, while still allowing the user to fully control the parameters
#   of the chosen road.
#
# REQUIRED SOURCES:
#   source("generate_KL_process.R")
#   source("generate_decomposition_process.R")
#
# Dependencies: base R (+ MASS indirectly through generate_decomposition_process)
# =============================================================================



#' Generate centred process through KL or decomposition road
#'
#' @param T_list        List of n numeric vectors. T_list((i)) contains the m_i
#'                      measurement times for subject i.
#' @param road          Character. Must be either "kl" or "decomposition".
#'
#' @param eigenfn_list  List of K vectorised eigenfunctions. Used only if
#'                      road = "kl".
#' @param eigenval_vec  Numeric vector of K positive eigenvalues. Used only if
#'                      road = "kl".
#' @param score_type    Character. Distribution family for KL scores.
#'                      Used only if road = "kl". Default: "gaussian".
#' @param score_params  Named list of parameters for score distribution.
#'                      Used only if road = "kl". Default: list().
#'
#' @param sigmaX_fn     Function. Standard deviation function sigma_X(t).
#'                      Used only if road = "decomposition".
#' @param corr_type     Character. Correlation family. Used only if
#'                      road = "decomposition".
#' @param corr_params   Named list of correlation parameters. Used only if
#'                      road = "decomposition". Default: list().
#' @param jitter        Small non-negative number added to the diagonal of
#'                      Sigma_i for numerical stability. Used only if
#'                      road = "decomposition". Default: 1e-10.
#'
#' @param seed          Integer or NULL. Optional seed for reproducibility.
#'
#' @return A list. Structure depends on the chosen road:
#'   - road = "kl":
#'       list(
#'         Z_list = ...,
#'         scores = ...
#'       )
#'   - road = "decomposition":
#'       list(
#'         Z_list = ...,
#'         Sigma_list = ...,
#'         config = ...
#'       )
#'
#' @examples
#' # ---------------------------------------------------------------------------
#' # Example 1: KL road
#' # ---------------------------------------------------------------------------
#' T_list <- list(c(0.1, 0.2, 0.3), c(0.55, 0.60, 0.72, 0.80))
#'
#' K <- 10
#' eigenfn_list <- lapply(1:K, function(k) {
#'   force(k)
#'   function(t) sqrt(2) * sin(2 * k * pi * t)
#' })
#' eigenval_vec <- 1 / (1:K)^2
#'
#' out_kl <- generate_centered_process(
#'   T_list = T_list,
#'   road = "kl",
#'   eigenfn_list = eigenfn_list,
#'   eigenval_vec = eigenval_vec,
#'   score_type = "gaussian",
#'   seed = 123
#' )
#'
#' # ---------------------------------------------------------------------------
#' # Example 2: Decomposition road
#' # ---------------------------------------------------------------------------
#' sigmaX_fn <- function(t) {
#'   sqrt( sqrt(t) * sqrt(1 - t) / 10 + 1 )
#' }
#'
#' out_dec <- generate_centered_process(
#'   T_list = T_list,
#'   road = "decomposition",
#'   sigmaX_fn = sigmaX_fn,
#'   corr_type = "matern",
#'   corr_params = list(theta1 = 0.5, theta2 = 1),
#'   seed = 123
#' )

#' @export
generate_centered_process <- function(T_list,
                                      road,
                                      eigenfn_list = NULL,
                                      eigenval_vec = NULL,
                                      score_type   = "gaussian",
                                      score_params = list(),
                                      sigmaX_fn    = NULL,
                                      corr_type    = NULL,
                                      corr_params  = list(),
                                      jitter       = 1e-10,
                                      seed         = NULL) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.list(T_list) || length(T_list) < 1)
    stop("T_list must be a non-empty list of numeric vectors.")
  
  if (missing(road))
    stop("You must specify road = 'kl' or road = 'decomposition'.")
  
  if (!is.character(road) || length(road) != 1 ||
      !road %in% c("kl", "decomposition")) {
    stop("road must be exactly one of: 'kl', 'decomposition'.")
  }
  
  # ---------------------------------------------------------------------------
  # ROAD = KL
  # ---------------------------------------------------------------------------
  if (road == "kl") {
    
    if (is.null(eigenfn_list) || is.null(eigenval_vec)) {
      stop("For road = 'kl', you must provide both eigenfn_list and eigenval_vec.")
    }
    
    out <- generate_KL_process(
      T_list        = T_list,
      eigenfn_list  = eigenfn_list,
      eigenval_vec  = eigenval_vec,
      score_type    = score_type,
      score_params  = score_params
    )
    
    return(out)
  }
  
  # ---------------------------------------------------------------------------
  # ROAD = DECOMPOSITION
  # ---------------------------------------------------------------------------
  if (road == "decomposition") {
    
    if (is.null(sigmaX_fn))
      stop("For road = 'decomposition', you must provide sigmaX_fn.")
    
    if (is.null(corr_type))
      stop("For road = 'decomposition', you must provide corr_type.")
    
    out <- generate_decomposition_process(
      T_list       = T_list,
      sigmaX_fn    = sigmaX_fn,
      corr_type    = corr_type,
      corr_params  = corr_params,
      seed         = seed,
      return_Sigma = TRUE,
      jitter       = jitter
    )
    
    return(out)
  }
}
