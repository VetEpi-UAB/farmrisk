#' Process and convert form version numbers
#'
#' @title Form Version Processing
#' @description Extracts and converts form version numbers from either character strings or data frames
#' into a numeric format. Handles version numbers in the format "vX_Y" where X is the major version
#' and Y is the subversion.
#'
#' @param data Either a character vector containing version strings or a data frame with a 'version' column
#' @return A numeric vector where each element is calculated as (major_version * 100 + subversion)
#' @export
#' @examples
#' # From character vector
#' form_version("v2_7")  # Returns 207
#' 
#' # From data frame
#' df <- data.frame(version = c("v1_0", "v2_3"))
#' form_version(df)  # Returns c(100, 203)

form_version <- function(data) {
  # Extract version text either from character input or data frame
  text <- if(is.character(data)) {
    data
  } else {
    data$version
  }
  
  # Initialize results vector
  result <- numeric(length(text))
  
  # Process each version string
  for(i in seq_along(text)) {
    text_i <- text[i]
    # Extract version numbers using regex
    matches <- regexec("v(\\d+)_(\\d+)", text_i)
    extracted <- regmatches(text_i, matches)[[1]]
    
    # Convert to numeric and calculate final version number
    version <- as.numeric(extracted[2])
    subversion <- as.numeric(extracted[3])
    result[i] <- version * 100 + subversion
  }
  
  return(result)
}

#' Handle backward compatibility for biosecurity survey versions
#'
#' @title Biosecurity Survey Backward Compatibility
#' @description Manages compatibility between different versions of bsg survey data,
#' applying necessary transformations to maintain consistency with the latest version.
#'
#' @param data A data frame containing bsg survey data
#' @return A modified data frame compatible with the latest version
#' @export
#' @examples
#' bsg <- read.csv("bsg.csv")
#' compatible_data <- bsg_backward(bsg)

bsg_backward <- function(data) {
  # Get latest version from forms directory
  version_list <- list.files("forms", pattern="^bsg_form.*\\.json$")
  last_version <- max(form_version(version_list))
  last_version_name <- gsub("bsg_form_|\\.json","",
                            version_list[form_version(version_list)==last_version])
  
  # Set backward compatibility threshold
  backward_version_name <- "v2_7_ciria"
  backward_version <- form_version(backward_version_name)
  
  # Get version of input data
  data_version <- form_version(data)
  
  # Check if data is already in latest version
  if(all(data_version == last_version)) {
    return(data)
  } else {
    message("Using old bsg version: ", 
            unique(data$version[!data_version==last_version]), 
            " (bsg last version: ", last_version_name, ")")
  }
  
  # Check compatibility
  if(any(data_version < backward_version)) {
    stop("bsg ", unique(data$version[data_version < backward_version]), 
         " not compatible. The current script requires mov_form over V2_3_ciria.")
  }
  
  # Version-specific transformations
  if("v2_7_ciria" %in% data$version) {
    # Fix outdoors/outdoor naming inconsistency
    data$outdoors_high_n[data$version %in% "v2_7_ciria"] <- 
      data$outdoor_high_n[data$version %in% "v2_7_ciria"]
    data$outdoors_low_n[data$version %in% "v2_7_ciria"] <- 
      data$outdoor_low_n[data$version %in% "v2_7_ciria"]
  }
  
  if(any(data_version < 214)) {
    # Update column names for consistency
    names(data) <- gsub("pasture_name", "pasture_id", names(data))
    names(data) <- gsub("veh_name", "veh_id", names(data))
  }
  
  return(data)
}

#' Handle backward compatibility for movement survey versions
#'
#' @title Movements Survey Backward Compatibility
#' @description Manages compatibility between different versions of mov survey data,
#' applying necessary transformations to maintain consistency with the latest version.
#'
#' @param data A data frame containing mov survey data
#' @return A modified data frame compatible with the latest version
#' @export
#' @examples
#' mov <- read.csv("mov.csv")
#' compatible_data <- mov_backward(mov)

mov_backward <- function(data) {
  # Get latest version from forms directory
  version_list <- list.files("forms", pattern="^mov_form.*\\.json$")
  last_version <- max(form_version(version_list))
  last_version_name <- gsub("mov_form_|\\.json","",
                            version_list[form_version(version_list)==last_version])
  
  # Set backward compatibility threshold
  backward_version_name <- "v2_3_ciria"
  backward_version <- form_version(backward_version_name)
  
  # Get version of input data
  data_version <- form_version(data)
  
  # Check if data is already in latest version
  if(all(data_version == last_version)) {
    return(data)
  } else {
    message("Using old mov version: ",
            unique(data$version[!data_version==last_version]),
            " (mov last version: ", last_version_name, ")")
  }
  
  # Check compatibility
  if(any(data_version < backward_version)) {
    stop("mov ", unique(data$version[data_version < backward_version]),
         " not compatible. The current script requires mov_form over V2_3_ciria.")
  }
  
  # Version-specific transformations
  if(any(data_version < 204)) {
    # Update column names
    names(data)[names(data) %in% "pasture_to_veh_mov_exclusive"] <- "pasture_to_veh_exclusive"
    names(data)[names(data) %in% "pasture_from_veh_mov_exclusive"] <- "pasture_from_veh_exclusive"
  }
  
  if(any(data_version < 205)) {
    # Fix exclusive flags
    data$pasture_to_veh_exclusive <- ifelse(data$pasture_to_veh_days_exclusive, "yes",
                                            data$pasture_to_veh_each_exclusive)
    data$pasture_from_veh_exclusive <- ifelse(data$pasture_from_veh_days_exclusive, "yes",
                                              data$pasture_from_veh_each_exclusive)
    
    # Update column names
    names(data) <- gsub("pasture_name", "pasture_id", names(data))
    names(data) <- gsub("veh_name", "veh_id", names(data))
  }
  
  if(any(data_version < 206)) {
    data["unk_plan" %in% data] <- "unk"
  }
  
  return(data)
}