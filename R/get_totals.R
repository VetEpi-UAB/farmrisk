#' @title Modular Monte Carlo Totals
#' @description Calculate probability totals, across events, levels or groups.
#' 
#' Calculate Combined Probability of Events
#' 
#' @description Calculates the probability of at least one event occurring from a list of events,
#' assuming independence between events.
#'
#' @param mcmodule mcmodule object containing the data and node list
#' @param mcnodes character vector of node names to combine
#' @param name optional custom name for the output node
#' @param prefix optional prefix for the output node name
#' @param data_name name/index of data to use if multiple datasets exist
#' @param summary logical, whether to calculate summary statistics
#' @param module_name optional module name override
#'
#' @return Updated mcmodule with new combined probability node
#' @export
#'
#' @examples
#' at_least_one(mcmodule = intro,
#'              mcnodes = c("no_purchase_inf_agg", "b_entry_agg"))
at_least_one <- function(mcmodule, mcnodes, name = NULL, prefix = NULL,
                         data_name = NULL, summary = calc_summary, module_name = NULL) {
  
  # Get appropriate data
  if (is.data.frame(mcmodule$data)) {
    data <- mcmodule$data
  } else {
    data_name <- ifelse(is.null(data_name), 1, data_name)
    data <- mcmodule$data[[data_name]]
  }
  
  if (is.null(module_name)) {
    module_name <- deparse(substitute(mcmodule))
  }
  
  p_all <- 0
  keys_names <- c()
  
  # Handle nodes with hierarchical grouping
  if (length(mcmodule$node_list[[mcnodes[1]]][["hg"]]) > 0) {
    for (i in 1:length(mcnodes)) {
      mc_name <- mcnodes[i]
      keys_type <- ifelse("agg_keys" %in% names(mcmodule$node_list[[mc_name]]),
                          "agg_keys", "keys")
      
      keys_names <- unique(c(keys_names, 
                             mcmodule$node_list[[mc_name]][[keys_type]]))
      p_i <- mcmodule$node_list[[mc_name]][["mcnode"]]
      
      # Match dimensions if needed
      if (dim(p_i)[3] < nrow(data)) {
        p_i <- mc_match(mcmodule, mc_name, data)
      }
      p_all <- 1 - ((1 - p_all) * (1 - p_i))
      
      hg_index <- data$hg
    }
    
  } else {
    # Handle nodes without hierarchical grouping
    if (!length(mcnodes) == 2) {
      stop("To aggregate mcnodes without hg index provide exactly two mc_nodes")
    }
    mc_name_x <- mcnodes[1]
    mc_name_y <- mcnodes[2]
    
    keys_type_x <- ifelse("agg_keys" %in% names(mcmodule$node_list[[mc_name_x]]),
                          "agg_keys", "keys")
    keys_type_y <- ifelse("agg_keys" %in% names(mcmodule$node_list[[mc_name_y]]),
                          "agg_keys", "keys")
    
    keys_names_x <- unique(c(keys_names, 
                             mcmodule$node_list[[mc_name_x]][[keys_type_x]]))
    keys_names_y <- unique(c(keys_names, 
                             mcmodule$node_list[[mc_name_y]][[keys_type_y]]))
    
    keys_names <- unique(intersect(keys_names_x, keys_names_y))
    
    p_xy <- mc_agg_match(mcmodule, mc_name_x, mc_name_y, keys_names)
    p_all <- 1 - ((1 - p_xy[[1]]) * (1 - p_xy[[2]]))
    data <- p_xy$index
    hg_index <- NULL
  }
  
  # Create new node name
  new_mc_name <- ifelse(is.null(name), mc_name, name)
  p_all_a_mc_name <- ifelse(is.null(name),
                            paste0(prefix, gsub("_tr$|_pi$", "", new_mc_name), "_all"),
                            new_mc_name)
  
  # Add node to module
  mcmodule$node_list[[p_all_a_mc_name]][["mcnode"]] <- p_all
  mcmodule$node_list[[p_all_a_mc_name]][["type"]] <- "total"
  mcmodule$node_list[[p_all_a_mc_name]][["param"]] <- c(mcnodes)
  mcmodule$node_list[[p_all_a_mc_name]][["inputs"]] <- c(mcnodes)
  mcmodule$node_list[[p_all_a_mc_name]][["description"]] <- 
    paste("Probability at least one of", mcnodes, "(assuming independence)")
  mcmodule$node_list[[p_all_a_mc_name]][["module"]] <- module_name
  mcmodule$node_list[[p_all_a_mc_name]][["keys"]] <- keys_names
  mcmodule$node_list[[p_all_a_mc_name]][["node_expression"]] <- 
    paste0("1-(", paste(paste("(1-", mcnodes, ")", sep = ""), collapse = "*"), ")")
  mcmodule$node_list[[p_all_a_mc_name]][["scenario"]] <- data$scenario_id
  mcmodule$node_list[[p_all_a_mc_name]][["hg"]] <- hg_index
  
  if (summary) {
    mcmodule$node_list[[p_all_a_mc_name]][["summary"]] <- 
      mc_summary(data = data,
                 mcnode = p_all,
                 mc_name = p_all_a_mc_name,
                 keys_names = keys_names)
  }
  
  return(mcmodule)
}


#' Calculate Total Across Levels
#' 
#' Calculates probabilities and expected counts at animal, farm and batch levels
#'
#' @param mcmodule mcmodule object containing the data
#' @param mcnodes character vector of node names to process
#' @param animals_n name of column with animal counts
#' @param farms_n name of column with farm counts
#' @param name optional custom name for output nodes
#' @param prefix optional prefix for output node names
#' @param all_mcnodes logical, whether to process all nodes
#' @param data_name name/index of data to use
#' @param summary logical, whether to calculate summaries
#'
#' @return Updated mcmodule with new probability nodes
#' @export
#'
#' @examples
#' get_totals(mcmodule = purchase_origin,
#'            mcnodes = c("a_inf", "a_inf_pi", "a_inf_tr"))
get_totals <- function(mcmodule, mcnodes, animals_n = "animals_n", 
                       farms_n = "farms_n", name = NULL, prefix = NULL,
                       all_mcnodes = TRUE, data_name = NULL, 
                       summary = calc_summary) {
  
  # Get appropriate data
  if (is.data.frame(mcmodule$data)) {
    data <- mcmodule$data
  } else {
    data_name <- ifelse(is.null(data_name), 1, data_name)
    data <- mcmodule$data[[data_name]]
  }
  module_name <- deparse(substitute(mcmodule))
  
  # Get animal counts as mcnodes
  if (animals_n %in% names(mcmodule$node_list)) {
    animals_mc <- mcmodule$node_list[[animals_n]][["mcnode"]]
    new <- FALSE
  } else if (animals_n %in% names(data)) {
    animals_mc <- mcdata(data[[animals_n]], type = "0", nvariates = nrow(data))
    new <- TRUE
  } else {
    create_mc_nodes(mcmodule$data, 
                    mctable = mcnode_admin[mcnode_admin$mcnode %in% animals_n,])
    assign("animals_mc", get(animals_n))
    new <- TRUE
  }
  
  # Add to node list
  mcmodule$node_list[[animals_n]][["mcnode"]] <- animals_mc
  mcmodule$node_list[[animals_n]][["type"]] <- "n_trials"
  mcmodule$node_list[[animals_n]][["module"]] <- module_name
  
  if (admin_wif) {
    if (new) {
      mcmodule$node_list[[animals_n]][["scenario"]] <- data$scenario_id
      mcmodule$node_list[[animals_n]][["hg"]] <- data$hg
    }
    
    if (dim(animals_mc)[3] < nrow(data)) {
      animals_mc <- mc_match(mcmodule, "animals_n", data)
    }
  }
  
  # Get farm counts as mcnodes
  if (farms_n %in% names(mcmodule$node_list)) {
    farms_mc <- mcmodule$node_list[[farms_n]][["mcnode"]]
    new <- FALSE
  } else if (farms_n %in% names(data)) {
    farms_mc <- mcdata(data[[farms_n]], type = "0", nvariates = nrow(data))
    new <- TRUE
  } else {
    create_mc_nodes(mcmodule$data,
                    mctable = mcnode_admin[mcnode_admin$mcnode %in% farms_n,])
    assign("farms_mc", get(farms_n))
    new <- TRUE
  }
  
  mcmodule$node_list[[farms_n]][["mcnode"]] <- farms_mc
  mcmodule$node_list[[farms_n]][["type"]] <- "n_trials"
  mcmodule$node_list[[farms_n]][["module"]] <- module_name
  
  if (admin_wif) {
    if (new) {
      mcmodule$node_list[[farms_n]][["scenario"]] <- data$scenario_id
      mcmodule$node_list[[farms_n]][["hg"]] <- data$hg
    }
    if (dim(farms_mc)[3] < nrow(data)) {
      farms_mc <- mc_match(mcmodule, "farms_n", data)
    }
  }
  
  prefix <- ifelse(is.null(prefix), "", paste0(prefix, "_"))
  
  # 1. Calculate combined probability for all nodes if requested
  if (all_mcnodes) {
    mcmodule <- at_least_one(mcmodule, mcnodes, name, prefix, data_name,
                             summary, module_name)
    
    new_mc_name <- ifelse(is.null(name),
                          paste0(prefix, mcnodes, "_all"), name)
    
    p_all_a_mc_name <- unique(gsub("_tr$|_pi$", "", new_mc_name))
    
    mcnodes <- c(mcnodes, p_all_a_mc_name)
  }
  
  # 2. Process each node
  for (i in 1:length(mcnodes)) {
    mc_name <- mcnodes[i]
    new_mc_name <- mc_name
    new_mc_name <- gsub(prefix, "", new_mc_name)
    keys_names <- mcmodule$node_list[[mc_name]][["keys"]]
    
    # Get animal-level probability
    p_a <- mcmodule$node_list[[mc_name]][["mcnode"]]
    
    if (dim(p_a)[3] < nrow(data)) {
      p_a <- mc_match(mcmodule, mc_name, data)
    }
    
    total_mcnodes <- c()
    
    # Calculate farm-level probability
    p_f <- 1 - (1 - p_a)^animals_mc
    p_f_mc_name <- paste0(prefix, gsub("^a", "f", new_mc_name))
    p_f_mc_name <- gsub("_a_", "_f_", p_f_mc_name)
    
    mcmodule$node_list[[p_f_mc_name]][["mcnode"]] <- p_f
    mcmodule$node_list[[p_f_mc_name]][["param"]] <- c(mc_name, animals_n)
    mcmodule$node_list[[p_f_mc_name]][["inputs"]] <- c(mc_name, animals_n)
    mcmodule$node_list[[p_f_mc_name]][["description"]] <- 
      paste("Probability of at least one", mc_name, "in a farm")
    mcmodule$node_list[[p_f_mc_name]][["node_expression"]] <- 
      paste0("1-(1-", mc_name, ")^", animals_n)
    
    total_mcnodes <- c(total_mcnodes, p_f_mc_name)
    
    # Calculate expected number at farm level
    n_f <- p_a * animals_mc
    n_f_mc_name <- paste0(prefix, gsub("^a", "f", new_mc_name), "_n")
    n_f_mc_name <- gsub("_a_", "_f_", n_f_mc_name)
    
    mcmodule$node_list[[n_f_mc_name]][["mcnode"]] <- n_f
    mcmodule$node_list[[n_f_mc_name]][["param"]] <- c(mc_name, animals_n)
    mcmodule$node_list[[n_f_mc_name]][["inputs"]] <- c(mc_name, animals_n)
    mcmodule$node_list[[n_f_mc_name]][["description"]] <- 
      paste("Expected number of", mc_name, "animals from a farm")
    mcmodule$node_list[[n_f_mc_name]][["node_expression"]] <- 
      paste0(mc_name, "*", animals_n)
    
    total_mcnodes <- c(total_mcnodes, n_f_mc_name)
    
    # Calculate batch-level probability
    p_b <- 1 - (1 - p_f)^farms_mc 
    p_b_mc_name <- paste0(prefix, gsub("^a", "b", new_mc_name))
    p_b_mc_name <- gsub("_a_", "_b_", p_b_mc_name)
    
    mcmodule$node_list[[p_b_mc_name]][["mcnode"]] <- p_b
    mcmodule$node_list[[p_b_mc_name]][["param"]] <- c(p_f_mc_name, farms_n)
    mcmodule$node_list[[p_b_mc_name]][["inputs"]] <- c(p_f_mc_name, farms_n)
    mcmodule$node_list[[p_b_mc_name]][["description"]] <- 
      paste("Probability of at least one", mc_name, "in the batch")
    mcmodule$node_list[[p_b_mc_name]][["node_expression"]] <- 
      paste0("1-(1-", p_f_mc_name, ")^", farms_n)
    total_mcnodes <- c(total_mcnodes, p_b_mc_name)
    
    # Calculate expected number at batch level
    n_b <- n_f * farms_mc
    n_b_mc_name <- paste0(prefix, gsub("^a", "b", new_mc_name), "_n")
    n_b_mc_name <- gsub("_a_", "_b_", n_b_mc_name)
    
    mcmodule$node_list[[n_b_mc_name]][["mcnode"]] <- n_b
    mcmodule$node_list[[n_b_mc_name]][["param"]] <- c(n_f_mc_name, farms_n)
    mcmodule$node_list[[n_b_mc_name]][["inputs"]] <- c(n_f_mc_name, farms_n)
    mcmodule$node_list[[n_b_mc_name]][["description"]] <- 
      paste("Expected number of", mc_name, "in the batch")
    mcmodule$node_list[[n_b_mc_name]][["node_expression"]] <- 
      paste0(n_f_mc_name, "*", farms_n)
    total_mcnodes <- c(total_mcnodes, n_b_mc_name)
    
    # Add metadata for all calculated nodes
    for (j in 1:length(total_mcnodes)) {
      total_mc_name <- total_mcnodes[j]
      total_mcnode <- mcmodule$node_list[[total_mc_name]][["mcnode"]]
      
      mcmodule$node_list[[total_mc_name]][["type"]] <- "total"
      mcmodule$node_list[[total_mc_name]][["module"]] <- module_name
      mcmodule$node_list[[total_mc_name]][["keys"]] <- keys_names
      mcmodule$node_list[[total_mc_name]][["scenario"]] <- data$scenario_id
      mcmodule$node_list[[total_mc_name]][["hg"]] <- data$hg
      
      if (summary) {
        mcmodule$node_list[[total_mc_name]][["summary"]] <- 
          mc_summary(data = data,
                     mcnode = total_mcnode,
                     mc_name = total_mc_name,
                     keys_names = keys_names)
      }
    }
  }
  
  mcmodule$modules <- unique(c(mcmodule$modules, module_name))
  return(mcmodule)
}

#' Calculate Aggregated Totals Across Groups
#'
#' Aggregates node values across specified grouping variables
#'
#' @param mcmodule mcmodule object containing the data
#' @param mcnode name of node to aggregate
#' @param keys_names grouping variables for aggregation
#' @param suffix suffix for output node name
#' @param data_name name/index of data to use
#' @param summary logical, whether to calculate summaries
#' @param agg_variates logical, whether to aggregate variates
#'
#' @return Updated mcmodule with new aggregated nodes
#' @export
#'
#' @examples
#' get_agg_totals(mcmodule = purchase_origin,
#'                mcnode = "a_inf_all")
get_agg_totals <- function(mcmodule, mcnode, 
                           keys_names = c("pathogen", "farm_id", "scenario_id"),
                           suffix = "agg", data_name = NULL,
                           summary = calc_summary,
                           agg_variates = FALSE) {
  
  # Get appropriate data
  if (is.data.frame(mcmodule$data)) {
    data <- mcmodule$data
  } else {
    data_name <- ifelse(is.null(data_name), 1, data_name)
    data <- mcmodule$data[[data_name]]
  }
  
  module_name <- deparse(substitute(mcmodule))
  mc_name <- mcnode
  mcnode <- mcmodule$node_list[[mc_name]][["mcnode"]]
  
  # Match dimensions if needed
  if (dim(mcnode)[3] < nrow(data)) {
    mcnode <- mc_match(mcmodule, mc_name, data)
  }
  
  agg_total_mc_name <- paste0(mc_name, "_", suffix)
  
  # Extract variates
  variates_list <- list()
  inv_variates_list <- list()
  for (i in 1:dim(mcnode)[3]) {
    variates_list[[i]] <- extractvar(mcnode, i)
    inv_variates_list[[i]] <- 1 - extractvar(mcnode, i)
  }
  
  # Create grouping index
  key_col <- data %>%
    select(all_of(keys_names)) %>%
    unite(everything(), col = "key", sep = ", ", remove = FALSE)
  
  key_levels <- unique(key_col$key)
  
  # Process each group
  for (i in 1:length(key_levels)) {
    index <- key_col$key %in% key_levels[i]
    
    if (grepl("_n$", mc_name)) {
      # Sum for counts
      total_lev <- Reduce("+", variates_list[index])
      mcmodule$node_list[[agg_total_mc_name]][["description"]] <- 
        paste("Number expected events by:", 
              paste(keys_names, collapse = ", "))
      mcmodule$node_list[[agg_total_mc_name]][["node_expression"]] <- 
        paste0(mc_name, "_1+", mc_name, "variate_2+... by:",
               paste(keys_names, collapse = ", "))
    } else {
      # Combine probabilities
      total_lev <- 1 - Reduce("*", inv_variates_list[index])
      mcmodule$node_list[[agg_total_mc_name]][["description"]] <- 
        paste("Probability at least one of the events happening by:",
              paste(keys_names, collapse = ", "))
      mcmodule$node_list[[agg_total_mc_name]][["node_expression"]] <- 
        paste0("1-((1-", mc_name, "_1)*(1-", mc_name, "_2)...) by:",
               paste(keys_names, collapse = ", "))
    }
    
    # Aggregate results
    if (agg_variates) {
      # One row per original variate
      agg_index <- mcdata(index, type = "0", nvariates = length(index))
      
      if (exists("total_agg")) {
        total_agg <- total_agg + agg_index * total_lev
      } else {
        total_agg <- agg_index * total_lev
      }
      
      new_keys_names <- mcmodule$node_list[[mc_name]][["keys"]]
      key_data <- data
    } else {
      # One row per result
      if (exists("total_agg")) {
        total_agg <- addvar(total_agg, total_lev)
      } else {
        total_agg <- total_lev
      }
      new_keys_names <- keys_names
      key_data <- unique(key_col)
    }
  }
  
  # Add aggregated node to module
  mcmodule$node_list[[agg_total_mc_name]][["mcnode"]] <- total_agg
  mcmodule$node_list[[agg_total_mc_name]][["type"]] <- "agg_total"
  mcmodule$node_list[[agg_total_mc_name]][["module"]] <- module_name
  mcmodule$node_list[[agg_total_mc_name]][["agg_data"]] <- key_levels
  mcmodule$node_list[[agg_total_mc_name]][["agg_keys"]] <- new_keys_names
  mcmodule$node_list[[agg_total_mc_name]][["keys"]] <- 
    mcmodule$node_list[[mc_name]][["keys"]]
  mcmodule$node_list[[agg_total_mc_name]][["inputs"]] <- mc_name
  
  if (summary) {
    mcmodule$node_list[[agg_total_mc_name]][["summary"]] <- 
      mc_summary(data = key_data,
                 mcnode = total_agg,
                 mc_name = agg_total_mc_name,
                 keys_names = new_keys_names)
  }
  
  mcmodule$modules <- unique(c(mcmodule$modules, module_name))
  return(mcmodule)
}