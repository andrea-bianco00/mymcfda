# =============================================================================
# build_noise_fn.R
# -----------------------------------------------------------------------------
# Builds a standardised noise function (mean 0, variance 1) based on the
# user's choice of distribution family. The returned function has signature
# noise_fn(m) -> numeric vector of length m.
#
# This function is a dispatcher: the user specifies a type and optional
# parameters, and receives back a ready-to-use noise_fn that can be
# passed directly to generate_noise().
#
# SUPPORTED FAMILIES:
#   "gaussian"  : standard normal N(0, 1). No parameters needed.
#   "t"         : Student's t, rescaled to variance 1.
#                 Parameters: df (degrees of freedom, must be > 2).
#                 Rescaling: t(df) / sqrt(df / (df - 2)).
#   "uniform"   : continuous uniform on [-sqrt(3), sqrt(3)], which has
#                 mean 0 and variance 1. No parameters needed.
#   "laplace"   : Laplace (double exponential), rescaled to variance 1.
#                 Generated as difference of two exponentials, rescaled.
#                 No parameters needed.
#
# The user can also pass a custom function directly via type = "custom",
# providing the function in params$fn. This allows any distribution
# the user wants, as long as it has mean 0 and variance 1.
#
# ROLE IN THE PIPELINE:
#   build_noise_fn  -->  noise_fn (ready-to-use function)
#   generate_noise(T_list, sigma0, noise_fn)  -->  noise_list
#
# INPUT:
#   type   : character. One of "gaussian", "t", "uniform", "laplace",
#            "custom".
#   params : named list of parameters for the chosen family.
#            "gaussian" : no parameters needed (default: list())
#            "t"        : list(df = ...), df > 2 required
#            "uniform"  : no parameters needed
#            "laplace"  : no parameters needed
#            "custom"   : list(fn = ...), where fn is a function with
#                         signature fn(m) -> numeric vector of length m,
#                         mean 0, variance 1.
#
# OUTPUT:
#   A function noise_fn(m) that generates m i.i.d. standardised random
#   variables with mean 0 and variance 1.
#
# EXAMPLES:
#   # Gaussian (default)
#   nfn <- build_noise_fn("gaussian")
#   nfn(5)  # 5 standard normal values
#
#   # Student's t with 5 degrees of freedom
#   nfn_t <- build_noise_fn("t", list(df = 5))
#   nfn_t(5)  # 5 values from rescaled t(5)
#
#   # Uniform
#   nfn_u <- build_noise_fn("uniform")
#   nfn_u(5)  # 5 values from Uniform(-sqrt(3), sqrt(3))
#
#   # Custom
#   my_fn <- function(m) sample(c(-1, 1), m, replace = TRUE)
#   nfn_c <- build_noise_fn("custom", list(fn = my_fn))
#
# Dependencies: base R only.
# =============================================================================


#' Build a standardised noise function
#'
#' @param type   Character. Distribution family.
#' @param params Named list of parameters for the chosen family.
#'
#' @return A function noise_fn(m) -> numeric vector of length m,
#'         mean 0, variance 1.

#' @export
build_noise_fn <- function(type = "gaussian", params = list()) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  valid_types <- c("gaussian", "t", "uniform", "laplace", "custom")
  
  if (!is.character(type) || length(type) != 1 || !(type %in% valid_types))
    stop("type must be one of: ", paste(valid_types, collapse = ", "), ".")
  
  if (!is.list(params))
    stop("params must be a named list.")
  
  # ---------------------------------------------------------------------------
  # BUILD noise_fn BASED ON TYPE
  # ---------------------------------------------------------------------------
  noise_fn <- switch(type,
                     
                     "gaussian" = {
                       function(m) rnorm(m)
                     },
                     
                     "t" = {
                       # extract df
                       df <- params$df
                       if (is.null(df))
                         stop("For type = 't', params must contain 'df' ",
                              "(degrees of freedom).")
                       if (!is.numeric(df) || length(df) != 1 || df <= 2)
                         stop("df must be a single number greater than 2 ",
                              "(required for finite variance).")
                       
                       # rescale factor: Var(t(df)) = df / (df - 2), so divide by sqrt of that
                       scale <- sqrt(df / (df - 2))
                       
                       function(m) rt(m, df = df) / scale
                     },
                     
                     "uniform" = {
                       # Uniform on [-sqrt(3), sqrt(3)] has mean 0 and variance 1
                       a <- sqrt(3)
                       function(m) runif(m, min = -a, max = a)
                     },
                     
                     "laplace" = {
                       # Laplace(0, 1/sqrt(2)) has mean 0 and variance 1
                       # Generated as difference of two Exp(1) scaled by 1/sqrt(2)
                       s <- 1 / sqrt(2)
                       function(m) s * (rexp(m, rate = 1) - rexp(m, rate = 1))
                     },
                     
                     "custom" = {
                       fn <- params$fn
                       if (is.null(fn) || !is.function(fn))
                         stop("For type = 'custom', params must contain 'fn' ",
                              "(a function with signature fn(m) -> numeric(m)).")
                       fn
                     }
  )
  
  # ---------------------------------------------------------------------------
  # VALIDATE THE BUILT FUNCTION
  # ---------------------------------------------------------------------------
  # Save RNG state BEFORE any random number generation
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  
  # Test that it returns the correct length
  test_out <- tryCatch(noise_fn(10),
                       error = function(e) {
                         stop("The built noise_fn failed on test call. ",
                              "Error: ", e$message)
                       })
  
  if (!is.numeric(test_out) || length(test_out) != 10)
    stop("The built noise_fn must return a numeric vector of the ",
         "requested length. Got length: ", length(test_out), ".")
  
  # Empirical check on mean and variance (large sample)
  set.seed(99999)
  big_sample <- noise_fn(10000)
  emp_mean <- mean(big_sample)
  emp_var  <- var(big_sample)
  
  # Restore RNG state
  if (!is.null(old_seed)) {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  } else {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  
  if (abs(emp_mean) > 0.1)
    warning("The noise function has empirical mean ", round(emp_mean, 4),
            " (expected ~0). Check that it is centred.")
  
  if (abs(emp_var - 1) > 0.2)
    warning("The noise function has empirical variance ", round(emp_var, 4),
            " (expected ~1). Check that it is standardised.")
  
  return(noise_fn)
}
