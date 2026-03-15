# =============================================================================
# prepare_data.R
# -----------------------------------------------------------------------------
# Validates and converts functional sparse data into the format required
# by all estimation functions in this project (T_list and Y_list).
#
# Accepted input formats:
#
#   Format 1 — long data frame:
#     A data frame with one row per observation, containing at least three
#     columns: subject id, observation time, and observed value.
#     Column names are specified by the user via id_col, t_col, y_col.
#
#   Format 2 — already lists:
#     T_list and Y_list are passed directly. The function only validates them.
#
# Output:
#   A list with two elements:
#     - T_list: list of length n, T_list[[i]] = vector of observation times
#     - Y_list: list of length n, Y_list[[i]] = vector of observed values
# =============================================================================


#' Validate and convert functional sparse data to T_list / Y_list format
#'
#' @param data    Either a data frame (long format) or NULL if T_list/Y_list
#'                are provided directly.
#' @param T_list  List of observation times. Used only if data is NULL.
#' @param Y_list  List of observed values. Used only if data is NULL.
#' @param id_col  Character string. Name of the subject id column in data.
#'                Default is "subject_id".
#' @param t_col   Character string. Name of the time column in data.
#'                Default is "t".
#' @param y_col   Character string. Name of the observed value column in data.
#'                Default is "Y".
#'
#' @return A list with components:
#'   - T_list: list of length n with observation times per subject.
#'   - Y_list: list of length n with observed values per subject.
#'
#' @examples
#' # Format 1: long data frame
#' df <- data.frame(
#'   subject_id = c(1,1,1, 2,2,2,2, 3,3),
#'   t          = c(0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 0.3, 0.4),
#'   Y          = c(1.2, 0.9, 1.1, 0.8, 1.0, 1.2, 0.9, 1.3, 1.1)
#' )
#' result <- prepare_data(data = df)
#'
#' # Format 2: already lists
#' T_list <- list(c(0.1, 0.2, 0.3), c(0.5, 0.6, 0.7, 0.8), c(0.3, 0.4))
#' Y_list <- list(c(1.2, 0.9, 1.1), c(0.8, 1.0, 1.2, 0.9), c(1.3, 1.1))
#' result <- prepare_data(T_list = T_list, Y_list = Y_list)

#' @export
prepare_data <- function(data   = NULL,
                         T_list = NULL,
                         Y_list = NULL,
                         id_col = "subject_id",
                         t_col  = "t",
                         y_col  = "Y") {
  
  # ---------------------------------------------------------------------------
  # FORMAT 1: long data frame
  # ---------------------------------------------------------------------------
  if (!is.null(data)) {
    
    # --- check data is a data frame -------------------------------------------
    if (!is.data.frame(data)) {
      stop("'data' must be a data frame.")
    }
    
    # --- check required columns exist -----------------------------------------
    missing_cols <- setdiff(c(id_col, t_col, y_col), names(data))
    if (length(missing_cols) > 0) {
      stop(sprintf(
        "Column(s) not found in data: %s.\n  Available columns: %s.",
        paste(missing_cols, collapse = ", "),
        paste(names(data), collapse = ", ")
      ))
    }
    
    # --- check column types ---------------------------------------------------
    if (!is.numeric(data[[t_col]])) {
      stop(sprintf("Column '%s' must be numeric.", t_col))
    }
    if (!is.numeric(data[[y_col]])) {
      stop(sprintf("Column '%s' must be numeric.", y_col))
    }
    
    # --- check for NA values --------------------------------------------------
    if (any(is.na(data[[id_col]]))) {
      stop(sprintf("Column '%s' contains NA values.", id_col))
    }
    if (any(is.na(data[[t_col]]))) {
      stop(sprintf("Column '%s' contains NA values.", t_col))
    }
    if (any(is.na(data[[y_col]]))) {
      stop(sprintf("Column '%s' contains NA values.", y_col))
    }
    
    # --- convert to T_list and Y_list -----------------------------------------
    subjects <- unique(data[[id_col]])
    n        <- length(subjects)
    
    T_list <- vector("list", n)
    Y_list <- vector("list", n)
    
    for (i in seq_len(n)) {
      idx        <- data[[id_col]] == subjects[i]
      T_list[[i]] <- data[[t_col]][idx]
      Y_list[[i]] <- data[[y_col]][idx]
    }
    
    cat(sprintf("Data converted from long format: %d subjects detected.\n", n))
    
    # ---------------------------------------------------------------------------
    # FORMAT 2: already lists
    # ---------------------------------------------------------------------------
  } else if (!is.null(T_list) && !is.null(Y_list)) {
    
    if (!is.list(T_list) || !is.list(Y_list)) {
      stop("T_list and Y_list must be lists.")
    }
    if (length(T_list) != length(Y_list)) {
      stop("T_list and Y_list must have the same length n.")
    }
    
    n <- length(T_list)
    cat(sprintf("List format detected: %d subjects.\n", n))
    
  } else {
    stop("Provide either 'data' (long data frame) or both 'T_list' and 'Y_list'.")
  }
  
  # ---------------------------------------------------------------------------
  # COMMON VALIDATION for both formats
  # ---------------------------------------------------------------------------
  m_vec <- sapply(T_list, length)
  
  # each subject must have at least 2 observations
  too_few <- which(m_vec < 2)
  if (length(too_few) > 0) {
    stop(sprintf(
      "Subject(s) %s have fewer than 2 observations. At least 2 are required.",
      paste(too_few, collapse = ", ")
    ))
  }
  
  # T_list and Y_list must have the same length within each subject
  inconsistent <- which(sapply(seq_len(n), function(i) {
    length(T_list[[i]]) != length(Y_list[[i]])
  }))
  if (length(inconsistent) > 0) {
    stop(sprintf(
      "Subject(s) %s have different lengths in T_list and Y_list.",
      paste(inconsistent, collapse = ", ")
    ))
  }
  
  # no NA values inside the lists
  na_T <- which(sapply(T_list, function(x) any(is.na(x))))
  na_Y <- which(sapply(Y_list, function(x) any(is.na(x))))
  if (length(na_T) > 0) {
    stop(sprintf("T_list contains NA values for subject(s): %s.",
                 paste(na_T, collapse = ", ")))
  }
  if (length(na_Y) > 0) {
    stop(sprintf("Y_list contains NA values for subject(s): %s.",
                 paste(na_Y, collapse = ", ")))
  }
  
  # observation times must be numeric
  non_numeric_T <- which(!sapply(T_list, is.numeric))
  if (length(non_numeric_T) > 0) {
    stop(sprintf("T_list contains non-numeric values for subject(s): %s.",
                 paste(non_numeric_T, collapse = ", ")))
  }
  
  # summary
  cat(sprintf("Observations per subject: min=%d, max=%d, mean=%.1f\n",
              min(m_vec), max(m_vec), mean(m_vec)))
  cat(sprintf("Time range: [%.4f, %.4f]\n",
              min(sapply(T_list, min)), max(sapply(T_list, max))))
  
  return(list(T_list = T_list, Y_list = Y_list))
}
