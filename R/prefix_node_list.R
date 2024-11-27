#' Add Prefix to Node Names
#' 
#' Adds a prefix to node_list names and all input nodes (non-prev)
#'
#' @param mcmodule An mcmodule or a node_list object
#' @param prefix String to add as prefix of the new mcmodule mcnodes (preserves previous prefixes)
#' @param rewrite_module Name of a module to rewrite (new name equals prefix) - used after copying modules
#'
#' @return Returns either:
#'   \itemize{
#'     \item A node_list if input was a node_list
#'     \item A model_object if input was an mcmodule
#'   }
#'
#' @examples
#' purchase <- add_prefix(
#'   mcmodule = transport_to,
#'   prefix = "transport_from",
#'   rewrite_module = "transport_to"
#' )
add_prefix <- function(mcmodule, prefix = NULL, rewrite_module = NULL) {
  # Extract node list
  node_list <- mcmodule$node_list
  
  # Get node names and modules
  node_names <- names(node_list)
  node_module <- sapply(node_list, "[[", "module")
  
  # Get inputs and their modules
  inputs <- unique(unlist(sapply(node_list, "[[", "inputs")))
  inputs_module <- unlist(sapply(node_list[inputs], "[[", "module"))
  
  # Set default prefix if none provided
  if (is.null(prefix)) {
    prefix <- deparse(substitute(mcmodule))
  }
  
  # Handle module rewriting if specified
  if (!is.null(rewrite_module)) {
    # Rename module
    node_module[rewrite_module %in% node_module] <- prefix
    inputs_module[rewrite_module %in% inputs_module] <- prefix
    
    # Remove prefix
    node_names <- gsub(paste0(rewrite_module, "_"), "", node_names)
    names(node_module) <- gsub(paste0(rewrite_module, "_"), "", names(node_module))
    names(inputs_module) <- gsub(paste0(rewrite_module, "_"), "", names(names(inputs_module)))
  }
  
  # Get unique modules
  modules <- unique(c(unlist(strsplit(names(unlist(mcmodule$model_expression)), split = "\\.")), prefix))
  
  # Add prefix to node names and inputs
  node_names[node_module[node_names] %in% modules] <- paste0(prefix, "_", node_names[node_module[node_names] %in% modules])
  new_inputs <- inputs
  new_inputs[inputs_module[inputs] %in% modules] <- paste0(prefix, "_", inputs[inputs_module[inputs] %in% modules])
  
  # Remove duplicated prefixes
  node_names <- gsub(paste0(prefix, "_", prefix), prefix, node_names)
  new_inputs <- gsub(paste0(prefix, "_", prefix), prefix, new_inputs)
  
  names(new_inputs) <- inputs
  
  # Update node list inputs and prefix
  for (i in 1:length(node_list)) {
    if (node_list[[i]][["module"]] %in% modules) {
      old_inputs <- node_list[[i]][["inputs"]]
      node_list[[i]][["inputs"]] <- new_inputs[old_inputs]
      node_list[[i]][["prefix"]] <- prefix
    }
  }
  
  names(node_list) <- node_names
  
  # Return appropriate object type
  if (class(mcmodule) == "mcmodule") {
    mcmodule$node_list <- node_list
    mcmodule$prefix <- prefix
    return(mcmodule)
  } else {
    return(node_list)
  }
}