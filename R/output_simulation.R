#' Save Monte Carlo Module Summary
#' 
#' Saves summary tables from Monte Carlo module nodes to CSV files
#' 
#' @param mcmodule List containing Monte Carlo module data
#' @param mcnodes Character vector of node names to process
#' @return None (saves files to disk)
save_mc_summary <- function(mcmodule, mcnodes) {
  # Get nodes with summaries
  summary_yes <- !sapply(sapply(mcmodule$node_list, "[[", "summary"), is.null)
  outputs <- mcnodes[mcnodes %in% names(summary_yes[summary_yes])]
  no_outputs <- mcnodes[!mcnodes %in% names(summary_yes[summary_yes])]
  
  # Warning for nodes without summaries
  if (length(no_outputs) > 0) {
    warning("No summary found in ", paste(no_outputs, collapse = ", "))
  }
  
  # Create output directory if needed
  output_folder_id <- file.path(output_path, farm_id)
  dir.create(output_folder_id, showWarnings = FALSE, recursive = TRUE)
  
  # Save summaries to CSV
  for (output in outputs) {
    file_path <- file.path(output_folder_id, paste0(output, ".csv"))
    write.table(mcmodule$node_list[[output]][["summary"]],
                file = file_path,
                row.names = FALSE,
                sep = ";")
  }
}


#' Save Table Summary
#'
#' Saves data tables to CSV files
#'
#' @param tables Named list of data frames to save
#' @return None (saves files to disk)
save_table_summary <- function(tables) {
  # Create output directory if needed
  output_folder_id <- file.path(output_path, farm_id)
  dir.create(output_folder_id, showWarnings = FALSE, recursive = TRUE)
  
  # Save tables to CSV using names from the list
  table_names <- names(tables)
  if (is.null(table_names)) {
    stop("Tables must be provided as a named list")
  }
  
  for (name in table_names) {
    file_path <- file.path(output_folder_id, paste0(name, ".csv"))
    write.table(tables[[name]],
                file = file_path,
                row.names = FALSE,
                sep = ";")
  }
}