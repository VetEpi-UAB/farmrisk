#' Get Monte Carlo Input Variables from CSV Files
#'
#' This function searches through CSV files in a specified directory to find variables
#' that are inputs for Monte Carlo nodes (mcnodes) defined in a Monte Carlo table.
#'
#' @param path Character string specifying the path to input files. Default is "input_files/".
#' @param mctable Data frame containing the list of Monte Carlo objects and their definitions.
#'   Default is mcnode_admin.
#'
#' @return Character vector of column names that are MC node inputs.
#' @export
#'
#' @examples
#' \dontrun{
#'   mc_inputs <- get_mc_inputs()
#'   mc_inputs <- get_mc_inputs(path = "my_files/", mctable = my_mc_table)
#' }
get_mc_inputs <- function(path = input_path, mctable = mcnode_admin) {
  # Load data files as list
  suppressMessages(
    data_list <- load_input_files(farm_id = farm_id,
                                  path = path, 
                                  create_df = FALSE)
    
  )
  
  # Get all column names from the data files
  col_names <- unlist(lapply(data_list, colnames))
  
  # Check if mc_inputs exist in data columns
  data_mc_inputs <- grepl(paste(paste0("\\<", mctable$mcnode, ".*"), 
                                collapse = "|"), 
                          col_names)
  
  # Return column names (of all files) that are mcnode inputs
  return(col_names[data_mc_inputs])
}
