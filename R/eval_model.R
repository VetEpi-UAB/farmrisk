#' Match Data or Monte Carlo Nodes By Scenario or Groups
#'
#' Functions for managing Monte Carlo nodes in a modular simulation framework.
#' 
#' Create Node List from Model Expression
#'
#' This function creates a list of nodes based on a given model expression, handling input,
#' output, and previous nodes with their respective properties and relationships.
#'
#' @param model_expression An R expression containing model calculations
#' @param param_names Optional named vector for parameter renaming
#' @param mctable Reference table for MC nodes, defaults to mcnode_admin
#' @param path Path to input files, defaults to "input_files/"
#'
#' @return A list of class "mcnode_list" containing node information
#' @export
#'
#' @examples
#' get_node_list(model_expression = indir_contact_expression)
get_node_list <- function(model_expression, param_names = NULL, 
                          mctable = mcnode_admin, path = "input_files/") {
  
  module <- gsub("_expression", "", deparse(substitute(model_expression)))
  
  # Initialize lists and vectors
  out_node_list <- list()
  all_nodes <- c()
  
  # Process output nodes from model expression
  for(i in 2:length(model_expression)) {
    node_name <- deparse(model_expression[[i]][[2]])
    node_expression <- paste0(deparse(model_expression[[i]][[3]]), collapse = "")
    out_node_list[[node_name]][["node_expression"]] <- node_expression
    
    # Extract input node names
    exp1 <- gsub("_", "975UNDERSCORE2023", node_expression)
    exp2 <- gsub("[^[:alnum:]]", ",", exp1)
    exp3 <- gsub("975UNDERSCORE2023", "_", exp2)
    exp4 <- c(strsplit(exp3, split = ",")[[1]])
    exp4 <- exp4[!exp4 %in% ""]
    exp5 <- suppressWarnings(exp4[is.na(as.numeric(exp4))])
    inputs <- unique(exp5)
    
    # Check NA removal
    na_rm <- any(inputs %in% "mcnode_na_rm")
    if(na_rm) {
      out_node_list[[node_name]][["na_rm"]] <- na_rm
    }
    
    # Filter function inputs
    inputs <- inputs[!inputs %in% "mcnode_na_rm"]
    fun_input <- c()
    if(length(inputs) > 0) {
      for(j in 1:length(inputs)) {
        input_name <- inputs[j]
        is_fun <- if(exists(input_name)) is.function(get(input_name)) else FALSE
        if(is_fun) fun_input <- c(fun_input, input_name)
      }
    }
    
    inputs <- inputs[!inputs %in% fun_input]
    all_nodes <- unique(c(all_nodes, node_name, inputs))
    
    # Set node type
    out_node_list[[node_name]][["type"]] <- 
      if(!grepl("[:alpha:]", node_expression)) "scalar" else "out_node"
    
    out_node_list[[node_name]][["inputs"]] <- inputs
    out_node_list[[node_name]][["module"]] <- module
  }
  
  # Rename parameters
  for(i in 1:length(all_nodes)) {
    all_nodes[i] <- if(all_nodes[i] %in% names(param_names)) {
      param_names[all_nodes[i]]
    } else {
      all_nodes[i]
    }
  }
  
  # Process input nodes
  in_node_list <- list()
  input_nodes <- all_nodes[all_nodes %in% as.character(mctable$mcnode)]
  
  if(!exists("all_inputs")) {
    all_inputs <<- get_mc_inputs(path)
  }
  
  if(length(input_nodes) > 0) {
    for(i in 1:length(input_nodes)) {
      node_name <- input_nodes[[i]]
      mc_row <- mctable[mctable$mcnode == node_name,]
      
      in_node_list[[node_name]][["type"]] <- "in_node"
      
      if(!is.na(mc_row$mc_func)) {
        in_node_list[[node_name]][["mc_func"]] <- as.character(mc_row$mc_func)
      }
      
      in_node_list[[node_name]][["description"]] <- as.character(mc_row$description)
      
      # Process input columns and files
      inputs_col <- all_inputs[grepl(paste0("\\<", node_name, "(\\>|_[^_]*\\>)"), 
                                     all_inputs)]
      if(length(inputs_col) > 0) {
        in_node_list[[node_name]][["inputs_col"]] <- inputs_col
      }
      
      input_file <- unique(gsub("[0-9]", "", names(inputs_col)))
      if(length(input_file) > 0) {
        in_node_list[[node_name]][["input_file"]] <- input_file
      }
      
      # Process keys
      keys_names <- unlist(strsplit(input_file, "_"))
      keys_map <- c(animal = "animal_category", 
                    region = "region_code",
                    visit = "visit_type",
                    point = "region_code",
                    bsg = "farm_id",
                    mov = "mov_id")
      
      for(key in names(keys_map)) {
        keys_names <- ifelse(keys_names == key, keys_map[key], keys_names)
      }
      
      if(admin_wif) {
        keys_names <- c(keys_names, "scenario_id")
      }
      
      if(length(keys_names) > 0) {
        in_node_list[[node_name]][["keys"]] <- keys_names
      }
      
      in_node_list[[node_name]][["module"]] <- module
      
      # Handle parameter renaming
      if(node_name %in% param_names) {
        node_name_expression <- names(param_names)[param_names %in% node_name]
        names(in_node_list)[names(in_node_list) %in% param_names] <- node_name_expression
        all_nodes[all_nodes %in% param_names] <- node_name_expression
      }
    }
  }
  
  # Process output node keys
  for(i in 1:length(out_node_list)) {
    inputs <- out_node_list[[i]][["inputs"]]
    if(length(inputs) > 0) {
      keys_names <- unique(unlist(lapply(inputs, function(x) {
        in_node_list[[x]][["keys"]]
      })))
      if(length(keys_names) > 0) {
        out_node_list[[i]][["keys"]] <- keys_names
      }
    }
  }
  
  # Process previous nodes
  prev_node_list <- list()
  prev_nodes <- all_nodes[!all_nodes %in% c(names(in_node_list), names(out_node_list))]
  
  if(length(prev_nodes) > 0) {
    for(i in 1:length(prev_nodes)) {
      node_name <- prev_nodes[i]
      is_fun <- if(exists(node_name)) is.function(get(node_name)) else FALSE
      if(!is_fun) {
        prev_node_list[[node_name]][["type"]] <- "prev_node"
      }
    }
  }
  
  # Combine all node lists
  node_list <- c(in_node_list, prev_node_list, out_node_list)
  class(node_list) <- "mcnode_list"
  
  return(node_list)
}

#' Get Nodes from Previous Module
#'
#' Retrieves nodes from a previous module and assigns them to the parent environment
#'
#' @param mcmodule An mcmodule or mcnode_list object
#' @param get_nodes Optional vector of node names to retrieve
#'
#' @return A subset of the node list containing requested nodes
#' @export
#'
#' @examples
#' get_previous_nodes(mcmodule)
get_previous_nodes <- function(mcmodule, get_nodes = NULL) {
  if(class(mcmodule) == "mcmodule") {
    node_list <- mcmodule$node_list
  } else if(class(mcmodule) == "mcnode_list") {
    node_list <- mcmodule
  } else {
    stop("mcmodule or mcnode_list object must be provided")
  }
  
  node_names <- names(node_list)
  node_names <- node_names[node_names %in% get_nodes]
  
  if(length(node_names) > 0) {
    for(i in 1:length(node_names)) {
      node_name <- node_names[i]
      assign(node_name, node_list[[node_name]][["mcnode"]], envir = parent.frame())
    }
  }
  
  return(node_list[node_names])
}

#' Evaluate Model Expression
#'
#' Evaluates a model expression and creates a comprehensive node list with results
#'
#' @param model_expression Model expression to evaluate
#' @param data Data frame containing input data
#' @param param_names Optional parameter renaming
#' @param prev_mcmodule Previous module reference
#' @param match_prev Match previous module dimensions
#' @param summary Calculate summary statistics
#' @param mctable Reference table for MC nodes
#' @param path Input files path
#' @param create_nodes Create new nodes flag
#' @param only_node_list Return only node list
#' @param analysis Perform model analysis
#'
#' @return An mcmodule object or node list
#' @export
#'
#' @examples
#' eval_model_expression(model_expression = origin_expression, 
#'                      data = purchase_origin_data)
eval_model_expression <- function(model_expression, data, param_names = NULL,
                                  prev_mcmodule = NULL, match_prev = FALSE,
                                  summary = calc_summary, mctable = mcnode_admin, 
                                  path = "input_files/", create_nodes = TRUE, 
                                  only_node_list = FALSE, analysis = model_analysis) {
  
  if(create_nodes & !is.null(data)) {
    create_mc_nodes(data)
  }
  
  # Process model expression
  if(is.list(model_expression)) {
    model_expression_list <- model_expression
  } else {
    expression_name <- gsub("_expression", "", deparse(substitute(model_expression)))
    model_expression_list <- list(model_expression)
    names(model_expression_list) <- expression_name
  }
  
  node_list <- list()
  modules <- c()
  mc_list <- list()
  
  # Process each expression
  for(i in 1:length(model_expression_list)) {
    model_expression_i <- model_expression_list[[i]]
    module <- names(model_expression_list)[[i]]
    
    node_list_i <- get_node_list(model_expression_i, param_names = param_names,
                                 mctable = mctable, path = path)
    
    prev_nodes <- names(node_list_i)[grepl("prev_node", node_list_i)]
    prev_nodes <- prev_nodes[!prev_nodes %in% names(node_list)]
    
    # Handle previous nodes
    if(length(prev_nodes) > 0) {
      if(is.null(prev_mcmodule)) {
        stop("prev_mcmodule for ", paste(prev_nodes, collapse = ", "), 
             " needed but not provided")
      } else {
        prev_mcmodule_list <- if(class(prev_mcmodule) == "mcmodule") {
          list(prev_mcmodule)
        } else {
          prev_mcmodule
        }
        
        for(j in 1:length(prev_mcmodule_list)) {
          prev_mcmodule_i <- prev_mcmodule_list[[j]]
          
          if(length(prev_nodes) > 0 & 
             any(!prev_nodes %in% names(prev_mcmodule_i$node_list))) {
            
            prefixes <- unlist(sapply(prev_mcmodule_i$node_list, "[[", "prefix"))
            new_names <- sapply(names(prefixes), function(x) {
              gsub(paste0("^", prefixes[x], "_"), "", x)
            })
            
            original_names <- names(prefixes)
            names(prefixes) <- new_names
            
            prev_nodes_names <- prev_nodes
            prev_nodes <- ifelse(prev_nodes %in% original_names,
                                 prev_nodes,
                                 paste0(prefixes[prev_nodes], "_", prev_nodes))
            names(prev_nodes) <- prev_nodes_names
            prev_param_names <- prev_nodes
          }
          
          prev_node_list_i <- get_previous_nodes(prev_mcmodule_i, 
                                                 get_nodes = prev_nodes)
          
          for(k in 1:length(prev_nodes)) {
            node_name <- prev_nodes[k]
            node_list_i[[node_name]] <- prev_node_list_i[[node_name]]
            
            if(admin_wif & match_prev) {
              match_prev_mcnode <- mc_match(prev_mcmodule, node_name, data)
              assign(node_name, match_prev_mcnode)
            }
          }
        }
      }
    }
    
    node_list <- c(node_list, node_list_i)
    
    new_param_names <- if(exists("prev_param_names")) {
      c(param_names, prev_param_names)
    } else {
      param_names
    }
    
    # Rename parameters
    if(!is.null(new_param_names)) {
      for(j in 1:length(new_param_names)) {
        expression_name <- names(new_param_names)[j]
        data_name <- new_param_names[j]
        
        if(exists(data_name)) {
          assign(expression_name, get(data_name))
        } else if(!is.null(prev_mcmodule$node_list[[data_name]])) {
          assign(expression_name, 
                 prev_mcmodule$node_list[[data_name]][["mcnode"]])
        }
      }
    }
    
    # Evaluate expression
    eval(model_expression_i)
    message("\n", module, " evaluated")
    
    mc_list_i <- list()
    
    # Complete node list
    for(j in 1:length(node_list)) {
      node_name <- names(node_list)[j]
      
      if(node_name %in% prev_nodes) next
      
      inputs <- node_list[[node_name]][["inputs"]]
      node_list[[node_name]][["param"]] <- inputs
      
      inputs[inputs %in% names(new_param_names)] <- 
        new_param_names[inputs[inputs %in% names(new_param_names)]]
      node_list[[node_name]][["inputs"]] <- inputs
      
      if(!is.null(prev_mcmodule) & 
         node_list[[node_name]][["type"]] == "out_node") {
        keys_names <- unique(unlist(lapply(inputs, function(x) {
          node_list[[x]][["keys"]]
        })))
        node_list[[node_name]][["keys"]] <- keys_names
      }
      
      mcnode <- get(node_name)
      
      if(!is.mcnode(mcnode) & is.numeric(mcnode)) {
        mcnode <- mcdata(mcnode, type = "0", nvariates = length(mcnode))
      }
      
      node_list[[node_name]][["mcnode"]] <- mcnode
      
      if(length(node_list[[node_name]][["module"]]) == 0 ||
         node_list[[node_name]][["module"]] %in% "model_i") {
        node_list[[node_name]][["module"]] <- module
      }
      
      modules <- unique(c(modules, node_list[[node_name]][["module"]]))
      
      if(summary & is.mc(mcnode)) {
        inputs_names <- node_list[[j]][["inputs"]]
        keys_names <- node_list[[j]][["keys"]]
        node_summary <- mc_summary(data = data, mcnode = mcnode,
                                   mc_name = node_name,
                                   keys_names = keys_names)
        node_list[[j]][["summary"]] <- node_summary
      }
      
      if(admin_wif) {
        node_list[[j]][["scenario"]] <- data$scenario_id
        node_list[[j]][["hg"]] <- data$hg
      }
      
      if(analysis) {
        for(k in 1:dim(mcnode)[3]) {
          mc_list_i_k_j <- mc(extractvar(mcnode, k), name = node_name)
          
          if(is.null(mc_list_i[k][[1]])) {
            mc_list_i[[k]] <- mc(mc_list_i_k_j)
          } else {
            mc_list_i[[k]] <- mc(mc_list_i[[k]], mc_list_i_k_j)
          }
        }
      }
    }
    
    mc_list[[module]] <- mc_list_i
  }
  
  # Clean up previous nodes
  node_list <- node_list[!(sapply(node_list, "[[", "type") == "prev_node")]
  
  if(only_node_list) {
    return(node_list)
  } else {
    mcmodule <- list(
      data = data,
      model_expression = model_expression,
      node_list = node_list,
      mc_list = mc_list,
      modules = modules
    )
    
    class(mcmodule) <- "mcmodule"
    
    message("\nmcmodule created (expressions: ",
            paste(names(model_expression), collapse = ", "), ")")
    
    return(mcmodule)
  }
}