# =============================================================================
# build_long_dataset.R
# -----------------------------------------------------------------------------
# Convert the simulated data produced by generate_full_snippet_data()
# into a long-format dataset containing only the observed data.
#
# OUTPUT STRUCTURE
#
#   id   time   Y
#
# This format is convenient for downstream estimation procedures.
# =============================================================================

build_long_dataset <- function(sim_data) {
  
  if (!is.list(sim_data)) {
    stop("Input must be the output of generate_full_snippet_data().")
  }
  
  if (!all(c("T_list", "Y_list") %in% names(sim_data))) {
    stop("sim_data must contain T_list and Y_list.")
  }
  
  T_list <- sim_data$T_list
  Y_list <- sim_data$Y_list
  
  n <- length(T_list)
  
  df_long <- do.call(
    rbind,
    lapply(seq_len(n), function(i) {
      
      data.frame(
        id   = i,
        time = T_list[[i]],
        Y    = Y_list[[i]]
      )
      
    })
  )
  
  rownames(df_long) <- NULL
  
  return(df_long)
}
