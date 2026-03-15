# =============================================================================
# eval_mean_at_times.R
# -----------------------------------------------------------------------------
# Evaluates a user-supplied mean function mu(t) at the measurement times
# T_i1, ..., T_im_i for each subject i.
#
# This function is a utility step between generate_snippet_data() (which
# produces T_list) and the process generation step (which needs mu(T_ij)
# to construct X_i(T_ij) = mu(T_ij) + sum_k xi_ik phi_k(T_ij)).
#
# ROLE IN THE PIPELINE:
#   generate_snippet_data  -->  T_list
#   eval_mean_at_times     -->  mu_values  (this file)
#   generate_KL_process    -->  process_values (centred, without mu)
#   generate_noise         -->  noise_values
#   assemble_observations  -->  Y_list = mu_values + process_values + noise
#
# The mean function mu is entirely user-specified: the user defines any
# R function mu_fn : [0,1] -> R and passes it here. This function does
# NOT impose any parametric form on mu.
#
# INPUT:
#   T_list : list of n numeric vectors. T_list[[i]] contains the m_i
#            measurement times for subject i, all in [0, 1].
#            Typically produced by generate_snippet_data()$T_list.
#   mu_fn  : function mapping a numeric vector t to a numeric vector of
#            the same length, representing mu(t). Must be vectorised,
#            i.e. mu_fn(c(0.1, 0.5, 0.9)) must return a vector of
#            length 3.
#
# OUTPUT:
#   A list of n numeric vectors. mu_values[[i]] is a vector of length m_i
#   containing mu(T_i1), ..., mu(T_im_i).
#
# EXAMPLES:
#   # Mean function from Lin & Wang (2022), Section 5
#   mu_fn <- function(t) 2 * t^2 * cos(2 * pi * t)
#
#   # Zero mean (centred process)
#   mu_fn <- function(t) rep(0, length(t))
#
#   # Exponential mean
#   mu_fn <- function(t) exp(t) / 2
#
# Dependencies: base R only.
# =============================================================================


#' Evaluate a mean function at measurement times
#'
#' @param T_list  List of n numeric vectors. T_list((i)) contains the m_i
#'                measurement times for subject i, all in (0, 1).
#' @param mu_fn   Function: numeric vector -> numeric vector of same length.
#'                The user-defined mean function mu(t). Must be vectorised.
#'
#' @return A list of n numeric vectors. mu_values((i)) contains
#'         mu(T_i1), ..., mu(T_im_i).
#'
#' @examples
#' # Setup
#' T_list <- list(c(0.1, 0.2, 0.3), c(0.5, 0.6), c(0.8, 0.9, 0.95, 0.99))
#' mu_fn  <- function(t) 2 * t^2 * cos(2 * pi * t)
#'
#' mu_vals <- eval_mean_at_times(T_list, mu_fn)
#' str(mu_vals)
#' # List of 3
#' #  $ : num (1:3) ...
#' #  $ : num (1:2) ...
#' #  $ : num (1:4) ...

#' @export
eval_mean_at_times <- function(T_list, mu_fn) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  if (!is.list(T_list) || length(T_list) < 1)
    stop("T_list must be a non-empty list of numeric vectors.")
  
  if (!is.function(mu_fn))
    stop("mu_fn must be a function.")
  
  # ---------------------------------------------------------------------------
  # CHECK VECTORISATION
  # ---------------------------------------------------------------------------
  # Test mu_fn on a small vector to verify it returns output of the
  # correct length. This catches common mistakes like using if/else
  # instead of ifelse(), which would return a scalar instead of a vector.
  test_input  <- c(0.0, 0.5, 1.0)
  test_output <- tryCatch(mu_fn(test_input),
                          error = function(e) {
                            stop("mu_fn failed when applied to a test vector ",
                                 "c(0, 0.5, 1). Error: ", e$message)
                          })
  
  if (!is.numeric(test_output))
    stop("mu_fn must return a numeric vector. Got: ", class(test_output)[1])
  
  if (length(test_output) != length(test_input))
    stop("mu_fn must be vectorised: mu_fn(t) must return a vector of the ",
         "same length as t. Input length: ", length(test_input),
         ", output length: ", length(test_output), ". ",
         "Hint: use ifelse() instead of if/else for conditional logic.")
  
  # ---------------------------------------------------------------------------
  # EVALUATE mu_fn AT EACH SUBJECT'S MEASUREMENT TIMES
  # ---------------------------------------------------------------------------
  n <- length(T_list)
  mu_values <- vector("list", n)
  
  for (i in seq_len(n)) {
    
    t_i <- T_list[[i]]
    
    # validate each element of T_list
    if (!is.numeric(t_i) || length(t_i) < 1)
      stop("T_list[[", i, "]] must be a non-empty numeric vector.")
    
    mu_values[[i]] <- mu_fn(t_i)
    
    # safety check: output length must match input length
    if (length(mu_values[[i]]) != length(t_i))
      stop("mu_fn returned ", length(mu_values[[i]]), " values for subject ", i,
           " but T_list[[", i, "]] has ", length(t_i), " time points. ",
           "mu_fn must be vectorised.")
  }
  
  return(mu_values)
}
