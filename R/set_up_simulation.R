#' Farm Risk Analysis Simulation Setup
#' 
#' @title Set up simulation for farm risk analysis
#' @description Initializes and configures a new simulation environment for farm risk analysis.
#' Sets up parameters, loads required input files, and processes biosecurity and movement data.
#' Assigns all simulation variables and states to the envir environment.
#'
#' @param user_id Character string. User identifier. Must be provided if farm_id is NULL.
#' @param farm_id Character. Farm identifier. Must be provided if user_id is NULL.
#' @param envir Environment. The environment where simulation variables will be assigned. Defaults to parent.frame().
#' @param ... Additional arguments passed to internal functions.
#'
#' @return None. All variables are assigned to envir environment
#'
#' @examples
#' # Basic usage
#' set_up_simulation("dev")
#' print(bsg)  # Access processed data directly
#' print(mov)  # Access movement data if it exists
#'
#' @export

set_up_simulation <- function(user_id=NULL,farm_id=NULL, envir = parent.frame(),...){
  
  # Validate that exactly one ID is provided
  if (is.null(user_id) && is.null(farm_id)) {
    stop("Must provide either user_id or farm_id")
  }
  if (!is.null(user_id) && !is.null(farm_id)) {
    stop("Provide only one ID: either user_id or farm_id, not both")
  }
  # Default settings
  defaults_set<-list(
    user_id=NULL,
    farm_id=NULL,
    risk_days=360,
    calc_summary = TRUE,
    admin_wif = TRUE,
    model_analysis = FALSE,
    alternative_params = NULL,
    alternative_formula = NULL,
    input_path=paste0(system.file("input_files/",package="farmrisk"),"/"),
    output_path=paste0(system.file("output_files/",package="farmrisk"),"/"),
    forms_path=paste0(system.file("forms/",package="farmrisk"),"/")
  )
  
  # Assign sefault settings
  for(i in 1:length(defaults_set)){
    if(!exists(names(defaults_set[i]))||is.null(get(names(defaults_set[i])))){
      assign(names(defaults_set[i]), defaults_set[[i]], envir = envir)
      assign(names(defaults_set[i]), defaults_set[[i]])
    }else{
      assign(names(defaults_set[i]), get(names(defaults_set[i])), envir = envir)
    }
  }
  
  # Display simulation settings
  message("\nsimulation_set: \n",
          " - calc_summary ", calc_summary,
          "\n - admin_wif: ", admin_wif,
          "\n - model_analysis: ", model_analysis,
          "\n - alternative_params: ", alternative_params,
          "\n - alternative_formula: ", deparse(alternative_formula),
          "\n - input_path: ", input_path,
          "\n - output_path: ", output_path,
          "\n - forms_path: ", forms_path
  )
  
  # Load and process input files
  load_input_files(user_id = user_id, farm_id = farm_id, envir = envir)
 
  # List of mc inputs
  assign("all_inputs",get_mc_inputs(),envir=envir)
  
  # Process biosecurity survey dates and assign to envir
  bsg <- process_bsg(envir$bsg)
  assign("bsg", bsg, envir = envir)
  
  # Process movement survey data if it exists
  if (exists("mov",envir=envir)) {
    mov <- process_mov(envir$mov, envir$risk_days)
    assign("mov", mov, envir = envir)
  }
  
  # No return needed as everything is in envir environment
  invisible(NULL)
}

#' Process biosecurity survey
#' @keywords internal
process_bsg <- function(bsg) {
  # Add current time if date is missing
  if (length(bsg$date) == 0) {
    bsg$date <- Sys.time()
  }
  
  # Convert dates to Unix timestamp
  if (!is.numeric(bsg$date)) {
    bsg$date <- as.numeric(as.POSIXct(bsg$date, 
                                      tryFormats = c("%d/%m/%Y","%Y-%m-%d")))
  }
  
  # Get latest survey
  bsg_index <- bsg$date == max(as.numeric(bsg$date))
  bsg_index <- cumsum(bsg_index) == max(cumsum(bsg_index))
  bsg <- bsg[bsg_index,]
  
  # Clean and update survey
  bsg$survey_id <- NULL
  bsg <- bsg_backward(bsg)
  bsg$scenario_id <- "0"
  
  message("bsg ", bsg$version, ": ", bsg$farm_id, " (", 
          as.POSIXct(bsg$date, origin = "1960-01-01"), ")")
  
  return(bsg)
}

#' Process movement survey data
#' @keywords internal
process_mov <- function(mov, risk_days) {
  # Handle dates
  if (length(mov$date) == 0) mov$date <- Sys.time()
  
  if (!is.numeric(mov$date)) {
    mov$date <- as.numeric(as.POSIXct(mov$date, 
                                      tryFormats = c("%d/%m/%Y","%Y-%m-%d")))
  }
  
  # Convert movement dates to days
  mov$purchase_single_date <- as.numeric(as.POSIXct(mov$purchase_single_date, 
                                                    tryFormats = c("%d/%m/%Y","%Y-%m-%d")))/86400
  mov$pasture_to_date <- as.numeric(as.POSIXct(mov$pasture_to_date,
                                               tryFormats = c("%d/%m/%Y","%Y-%m-%d")))/86400
  mov$pasture_from_date <- as.numeric(as.POSIXct(mov$pasture_from_date, 
                                                 tryFormats = c("%d/%m/%Y","%Y-%m-%d")))/86400
  
  # Filter by risk days
  last_date <- max(c(mov$purchase_single_date, mov$pasture_to_date, 
                     mov$pasture_from_date), na.rm = TRUE)
  first_date <- last_date - risk_days
  
  mov_filter <- (mov$purchase_single_date > first_date) | 
    (mov$pasture_to_date >= first_date)
  
  mov$pasture_from_date <- ifelse(mov$pasture_from_date < first_date, 
                                  first_date, mov$pasture_from_date)
  
  mov_filter <- ifelse(is.na(mov_filter), FALSE, mov_filter)
  
  # Display filtered movements
  if (any(!mov_filter)) {
    message("Some movements are prior to the risk days to calculate (", risk_days,")")
    message("The following movements are not included in the simulation: ",
            paste(mov[!mov_filter,c("mov_id")], collapse = ", "))
  }
  
  # Update movement data
  mov <- mov[mov_filter,]
  mov <- mov[!duplicated(mov$mov_id, fromLast = TRUE),]
  mov <- mov_backward(mov)
  mov$survey_id <- NULL
  mov$scenario_id <- "0"
  
  message("mov ", paste(unique(mov$version), collapse=", "), " (",
          nrow(mov), "): ", paste(mov$mov_id, collapse=", "))
  
  return(mov)
  
}