#' Load and Process Input Files for Farm Risk Analysis
#'
#' @title Load input files for farm risk analysis simulation
#' @description Loads and processes CSV and JSON files containing farm risk analysis data.
#' Files can be loaded from either admin or user-specific directories. Handles both
#' biosecurity (bsg) and movement (mov) data files.
#'
#' @param user_id Character string. Optional user identifier to load specific user files.
#' @param path Character string. Path to input files directory. Default: "input_files/"
#' @param from_json Logical. Whether to process JSON files if CSVs not found. Default: TRUE
#' @param create_df Logical. Whether to create data frames from files. Default: TRUE
#'
#' @return List containing loaded data frames
#'
#' @details
#' The function performs the following steps:
#' 1. Validates input path
#' 2. Loads admin files from admin directory
#' 3. If user_id provided:
#'    - Loads user-specific files
#'    - Checks for file modifications
#'    - Converts JSON to CSV if needed
#' 4. Processes all files and creates data frames
#'
#' @examples
#' load_input_files(user_id = "user123")
#' print(bsg)  # Access processed data directly
#'
#' @importFrom utils read.csv
#' @export

load_input_files <- function(user_id = NULL, 
                             farm_id = NULL,
                             path = input_path, 
                             from_json = TRUE, 
                             create_df = TRUE,
                             envir=parent.frame()) {
  
  # Validate path
  if (!file.exists(path)) stop("Path does not exist")
  
  # Process admin files
  admin_path <- paste0(path, "admin/")
  admin_names <- list.files(path = admin_path, pattern = "*.csv", recursive = TRUE)
  admin_files <- paste(admin_path, admin_names, sep = "")
  
  if (!is.null(user_id)) {
    # Read user configuration
    user_config <- read.csv(
      paste0(input_path,"user/id/", user_id, ".csv"),
      sep = ";",
      stringsAsFactors = FALSE,
      header = FALSE
    )
    
    # Store user information in envir environment
    assign("farm_id", user_config[1,1], envir = envir)
    assign("risk_days", user_config[1,2] * 30, envir = envir)

    
    # Process user files
    farm_id<-user_config[1,1]
    user_path <- paste0(path, "user/", farm_id, "/")
    user_files_info <- file.info(paste0(user_path, list.files(user_path)))
    
    # Check for file modifications
    process_json <- check_file_modifications(user_path, user_files_info)
    
    if (process_json && from_json) {
      process_json_files(user_path, envir)
    }
  }
    
    

  if(!is.null(farm_id)){
    #Stop if farm data was requested and not found
    if(!farm_id%in%list.files(path = paste0(path, "user/"))) stop(farm_id," not found in ",paste0(path, "user/"))
                                                                 
    # Process admin files
    user_path <- paste0(path, "user/",farm_id,"/")
    user_names <- list.files(path = user_path, pattern = "*.csv", recursive = TRUE)
    user_files <- paste(user_path, user_names, sep = "")
    
    file_names <- c(admin_names, user_names)
    file_list <- c(admin_files, user_files)
  }else{
    file_names <- admin_names
    file_list <- admin_files
  }
    
  
  # Load all CSV files
  data_list <- lapply(
    file_list[grepl(".csv$", file_list)],
    read.csv,
    sep = ";",
    stringsAsFactors = TRUE,
    na.strings = c(NA, "NA", "")
  )
  
  names(data_list) <- sub(".csv", "", file_names)
  
  # Assign data to envir environment
  if (length(data_list) > 0) {
    for (name in names(data_list)) {
      assign(name, data_list[[name]], envir = envir)
    }
    message("\nAdmin data loaded")
    if(is.null(user_id)){
      message("\nUser data loaded: ", farm_id)
    }else{
      message("\nUser data loaded: ", farm_id," (user: ",user_id,")")
    }
    
  }
  
  invisible(data_list)
}

# Helper functions for load_input_files

#' Check if files need to be processed from JSON
#' 
#' @param user_path Character string. Path to user directory
#' @param user_files_info File information from file.info()
#' @return Logical. TRUE if JSON processing is needed
check_file_modifications <- function(user_path, user_files_info) {
  # Check bsg file modifications
  bsg_change <- FALSE
  if (file.exists(paste0(user_path, "bsg.csv"))) {
    bsg_change <- user_files_info[paste0(user_path, "bsg.csv"), "mtime"] <
      max(user_files_info[grepl("bsg.json$", row.names(user_files_info)), "mtime"])
  }
  
  # Check mov file modifications
  mov_change <- FALSE
  if (!is.na(user_files_info[paste0(user_path, "mov.csv"), "mtime"])) {
    mov_change <- user_files_info[paste0(user_path, "mov.csv"), "mtime"] <
      max(user_files_info[grepl("mov_.*json$", row.names(user_files_info)), "mtime"])
  }
  
  return(bsg_change || mov_change)
}

#' Process JSON files to CSV
#' 
#' @param user_path Character string. Path to user directory
#' @param envir Environment. envir environment to store data
#' @importFrom jsonlite fromJSON
process_json_files <- function(user_path, envir) {
  # Process bsg file
  bsg_path <- paste0(user_path, "bsg.json")
  if (file.exists(bsg_path)) {
    bsg <- answer_json_to_df(bsg_path)
    assign("bsg", bsg, envir = envir)
    write.csv(bsg, 
              paste0(user_path, "bsg.csv"), 
              row.names = FALSE, 
              sep = ";")
    message("bsg.csv updated")
  }
  
  # Process mov files
  mov_files <- list.files(user_path, pattern = "mov_.*json$")
  if (length(mov_files) > 0) {
    mov <- do.call(rbind, lapply(mov_files, function(f) {
      answer_json_to_df(paste0(user_path, f))
    }))
    assign("mov", mov, envir = envir)
    write.csv(mov, 
              paste0(user_path, "mov.csv"), 
              row.names = FALSE, 
              sep = ";")
    message("mov.csv updated")
  }
}