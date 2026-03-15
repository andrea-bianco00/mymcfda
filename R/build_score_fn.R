# =============================================================================
# build_score_fn.R
# -----------------------------------------------------------------------------
# Builds a standardised score function (mean 0, variance 1) based on the
# user's choice of distribution family. The returned function has signature
# score_fn(m) -> numeric vector of length m.
#
# This function is a dispatcher: the user specifies a type and optional
# parameters, and receives back a ready-to-use score_fn. The scores
# xi_{ik} in the Karhunen-Loève expansion are then generated as:
#
#   xi_{ik} = sqrt(lambda_k) * score_fn(1)
#
# so that E[xi_{ik}] = 0 and Var(xi_{ik}) = lambda_k.
#
# By default, score_fn generates standard normal variables, giving a
# Gaussian process (as in Lin & Wang 2022, Section 5). The user can
# choose other distributions to test robustness of the estimators
# against non-Gaussianity.
#
# SUPPORTED FAMILIES:
#   "gaussian"  : standard normal N(0, 1). No parameters needed.
#                 This is the default, as in Lin & Wang (2022).
#   "t"         : Student's t, rescaled to variance 1.
#                 Parameters: df (degrees of freedom, must be > 2).
#   "uniform"   : continuous uniform on [-sqrt(3), sqrt(3)], which has
#                 mean 0 and variance 1. No parameters needed.
#   "laplace"   : Laplace (double exponential), rescaled to variance 1.
#                 No parameters needed.
#   "custom"    : user-supplied function via params$fn.
#
# ROLE IN THE PIPELINE:
#   build_score_fn      -->  score_fn (ready-to-use function)
#   generate_KL_process uses score_fn to generate xi_{ik}
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
#   A function score_fn(m) that generates m i.i.d. standardised random
#   variables with mean 0 and variance 1.
#
# EXAMPLES:
#   # Gaussian (default, as in Lin & Wang 2022)
#   sfn <- build_score_fn("gaussian")
#   sfn(5)  # 5 standard normal values
#
#   # Student's t with 5 degrees of freedom
#   sfn_t <- build_score_fn("t", list(df = 5))
#   sfn_t(5)  # 5 values from rescaled t(5)
#
#   # Uniform
#   sfn_u <- build_score_fn("uniform")
#   sfn_u(5)  # 5 values from Uniform(-sqrt(3), sqrt(3))
#
# Dependencies: base R only.
# =============================================================================


#' Build a standardised score function for KL expansion
#'
#' @param type   Character. Distribution family.
#' @param params Named list of parameters for the chosen family.
#'
#' @return A function score_fn(m) -> numeric vector of length m,
#'         mean 0, variance 1.

build_score_fn <- function(type = "gaussian", params = list()) {
  
  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  valid_types <- c("gaussian", "t", "uniform", "laplace", "custom")
  
  if (!is.character(type) || length(type) != 1 || !(type %in% valid_types))
    stop("type must be one of: ", paste(valid_types, collapse = ", "), ".")
  
  if (!is.list(params))
    stop("params must be a named list.")
  
  # ---------------------------------------------------------------------------
  # BUILD score_fn BASED ON TYPE
  # ---------------------------------------------------------------------------
  score_fn <- switch(type,
                     
                     "gaussian" = {
                       function(m) rnorm(m)
                     },
                     
                     "t" = {
                       df <- params$df
                       if (is.null(df))
                         stop("For type = 't', params must contain 'df' ",
                              "(degrees of freedom).")
                       if (!is.numeric(df) || length(df) != 1 || df <= 2)
                         stop("df must be a single number greater than 2 ",
                              "(required for finite variance).")
                       
                       scale <- sqrt(df / (df - 2))
                       function(m) rt(m, df = df) / scale
                     },
                     
                     "uniform" = {
                       a <- sqrt(3)
                       function(m) runif(m, min = -a, max = a)
                     },
                     
                     "laplace" = {
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
  
  test_out <- tryCatch(score_fn(10),
                       error = function(e) {
                         stop("The built score_fn failed on test call. ",
                              "Error: ", e$message)
                       })
  
  if (!is.numeric(test_out) || length(test_out) != 10)
    stop("The built score_fn must return a numeric vector of the ",
         "requested length. Got length: ", length(test_out), ".")
  
  # Empirical check on mean and variance (large sample)
  set.seed(99998)
  big_sample <- score_fn(10000)
  emp_mean <- mean(big_sample)
  emp_var  <- var(big_sample)
  
  # Restore RNG state
  if (!is.null(old_seed)) {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  } else {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  
  if (abs(emp_mean) > 0.1)
    warning("The score function has empirical mean ", round(emp_mean, 4),
            " (expected ~0). Check that it is centred.")
  
  if (abs(emp_var - 1) > 0.2)
    warning("The score function has empirical variance ", round(emp_var, 4),
            " (expected ~1). Check that it is standardised.")
  
  return(score_fn)
}
