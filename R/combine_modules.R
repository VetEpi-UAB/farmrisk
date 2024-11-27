#' Combine Two Modules
#'
#' This function combines two mcmodules into a single mcmodule, merging their data,
#' model expressions, node lists, and other components.
#'
#' @param mcmodule_x First module to combine
#' @param mcmodule_y Second module to combine
#' 
#' @return A combined mcmodule object containing elements from both input modules
#' 
#' @examples
#' combine_modules(mcmodule_x = purchase, mcmodule_y = fattening)
#' 
#' @export
combine_modules <- function(mcmodule_x, mcmodule_y) {
  mcmodule <- list()
  
  # Extract names of input modules
  name_x <- deparse(substitute(mcmodule_x))
  name_y <- deparse(substitute(mcmodule_y))
  
  # Check data structure of first module
  if (!is.data.frame(mcmodule_x$data)) {
    data_names_x <- names(mcmodule_x$data)
    islist_x <- TRUE
  } else {
    data_names_x <- name_x
    islist_x <- FALSE
  }
  
  # Check data structure of second module
  if (!is.data.frame(mcmodule_y$data)) {
    data_names_y <- names(mcmodule_y$data)
    islist_y <- TRUE
  } else {
    data_names_y <- name_y
    islist_y <- FALSE
  }
  
  # Combine data based on structure
  if (identical(mcmodule_x$data, mcmodule_y$data)) {
    mcmodule$data <- mcmodule_x$data
  } else if (!any(islist_x, islist_y)) {
    mcmodule$data <- list(mcmodule_x$data, mcmodule_y$data)
    names(mcmodule$data) <- c(data_names_x, data_names_y)
  } else if (all(islist_x, islist_y)) {
    mcmodule$data <- c(mcmodule_x$data, mcmodule_y$data)
    names(mcmodule$data) <- c(data_names_x, data_names_y)
  } else if (islist_x) {
    mcmodule$data <- c(mcmodule_x$data, list(data_name_y = mcmodule_y$data))
    names(mcmodule$data) <- c(data_names_x, data_names_y)
  } else if (islist_y) {
    mcmodule$data <- c(list(data_name_x = mcmodule_x$data), mcmodule_y$data)
    names(mcmodule$data) <- c(data_names_x, data_names_y)
  } else {
    stop("An error occurred combining data:", paste(data_names_x), " and ", paste(data_names_y))
  }
  
  # Combine model expressions
  mcmodule$model_expression <- list(
    mcmodule_x$model_expression,
    mcmodule_y$model_expression
  )
  names(mcmodule$model_expression) <- c(name_x, name_y)
  
  # Combine node lists and modules
  mcmodule$node_list <- c(mcmodule_x$node_list, mcmodule_y$node_list)
  mcmodule$modules <- unique(c(mcmodule_x$modules, mcmodule_y$modules))
  
  # Combine mc_list if they exist
  if (!is.null(mcmodule_x$mc_list)) {
    mcmodule$mc_list <- list(
      mcmodule_x$mc_list,
      mcmodule_y$mc_list
    )
    names(mcmodule$mc_list) <- c(name_x, name_y)
  }
  
  # Set class and return
  class(mcmodule) <- "mcmodule"
  return(mcmodule)
}