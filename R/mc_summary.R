#' @title Monte Carlo Node Keys and Summaries 
#' @description Collection of functions for handling Monte Carlo Nodes using their metadata found in their mcmodule
#' @export

#' Get mcnode admin inputs key names
#' @param x mcnode or character string
#' @return character vector of key names
#' @export
#' @examples
#' input_keys(mcnode = inf_dc)
input_keys <- function(x) {
  if (!exists("all_inputs")) {
    all_inputs <<- get_mc_inputs()
  }
  
  x_name <- if(is.character(x)) x else deparse(substitute(x))
  mc_inputs <- all_inputs[grepl(paste0("\\<", x_name, "(\\>|_[^_]*\\>)"), all_inputs)]
  
  if (length(mc_inputs) > 0) {
    input_df <- unique(gsub("[0-9]", "", names(mc_inputs)))
    input_col <- TRUE
  } else {
    input_df <- unique(gsub("[0-9]", "", names(all_inputs)))
    if (!x_name %in% input_df) {
      message(paste0(x_name, " is not an input file, keys not found."))
      return(NA)
    }
    input_df <- x_name
    input_col <- FALSE
  }
  
  # Key name mappings
  key_mappings <- c(
    animal = "animal_category",
    region = "region_code", 
    visit = "visit_type",
    point = "contact_point",
    density = "density_level",
    bsg = "farm_id",
    mov = "mov_id"
  )
  
  keys_names <- unlist(strsplit(input_df, "_"))
  for (old_key in names(key_mappings)) {
    keys_names[keys_names == old_key] <- key_mappings[old_key]
  }
  
  if (admin_wif & input_col) {
    keys_names <- c(keys_names, "scenario_id")
  }
  
  return(keys_names)
}

#' Get mcnode summary statistics
#' @param data data frame containing mnode inputs
#' @param mcnode mcnode object
#' @param sep_keys logical, separate keys in different columns
#' @param na_value value to exclude from summary
#' @param keys_names vector of key names for grouping
#' @param mc_name character, name of mcnode
#' @param by_keys logical, whether to summarize by keys
#' @return summary data frame
#' @export
#' @examples
#' mc_summary(mcnode = h_prev)
mc_summary <- function(data, mcnode, sep_keys = TRUE, na_value = NA, 
                       keys_names = NULL, mc_name = NULL, by_keys = FALSE) {
  
  if (!is.null(keys_names) && !all(keys_names %in% names(data))) {
    stop("keys_names (", paste(keys_names[!keys_names %in% names(data)], 
                               collapse = ", "), ") must appear in data column names")
  }
  
  if (!is.mcnode(mcnode)) {
    stop(mcnode, " must be a mc object")
  }
  
  mc_name <- if(is.null(mc_name)) deparse(substitute(mcnode)) else mc_name
  keys_names <- if(is.null(keys_names)) input_keys(mc_name) else keys_names
  
  keys <- if(length(keys_names) > 0 && any(keys_names %in% names(data))) {
    data[names(data) %in% keys_names]
  } else {
    data.frame(variate = 1:nrow(data))
  }
  
  if (!sep_keys) {
    keys <- unite(keys, col = keys, 1:ncol(keys), sep = ", ")
    keys_groups <- c("mc_name", "keys")
  } else {
    keys_groups <- c("mc_name", names(keys))
  }
  
  summary_l <- summary(mcnode)[[1]]
  if (!is.list(summary_l)) summary_l <- list(summary_l)
  
  summary_names <- colnames(summary_l[[1]])
  summary_df <- data.frame(matrix(unlist(summary_l), nrow = length(summary_l), 
                                  byrow = TRUE))
  names(summary_df) <- summary_names
  summary_df <- cbind(mc_name, keys, summary_df)
  
  if (by_keys) {
    filter_na_value <- if(is.na(na_value)) {
      if ("mean" %in% names(summary_df)) !is.na(summary_df$mean) 
      else !is.na(summary_df$NoUnc)
    } else {
      if ("mean" %in% names(summary_df)) !summary_df$mean == na_value 
      else !summary_df$NoUnc == na_value
    }
    
    summary_df <- summary_df[filter_na_value,] %>%
      group_by(across(all_of(keys_groups))) %>%
      summarise_all(mean) %>%
      mutate(across(where(is.numeric), my_round, digits = 4))
  }
  
  return(summary_df)
}

#' Get mcnode summary keys
#' @param mcsummary data frame from mc_summary()
#' @return vector of key names
#' @export
#' @examples
#' mc_summary_keys(mcsummary(inf_dc))
mc_summary_keys <- function(mcsummary) {
  if ("mean" %in% names(mcsummary)) {
    names(mcsummary)[2:(match("mean", names(mcsummary)) - 1)]
  } else {
    names(mcsummary)[2:(match("NoUnc", names(mcsummary)) - 1)]
  }
}

#' Include summary and keys in node_list
#' @param mcmodule mc module object
#' @param data data frame with mc inputs
#' @param node_list list of nodes
#' @return updated node_list
#' @export
#' @examples
#' node_list_summary(data = purchase_origin_data, node_list = node_list)
node_list_summary <- function(mcmodule = NULL, data = NULL, node_list = NULL) {
  if (!is.null(mcmodule)) {
    data <- mcmodule$data
    node_list <- mcmodule$node_list
  } else if (is.null(data) || is.null(node_list)) {
    stop("data containing mc_inputs and node_list must be provided")
  }
  
  for (i in seq_along(node_list)) {
    node_name <- names(node_list)[i]
    inputs_names <- node_list[[i]][["inputs"]]
    mcnode <- node_list[[i]][["mcnode"]]
    
    if (is.null(inputs_names)) {
      node_summary <- mc_summary(data = data, mcnode = mcnode, 
                                 mc_name = node_name)
    } else {
      keys_names <- unique(unlist(lapply(inputs_names, function(x) {
        node_list[[x]][["keys"]]
      })))
      
      node_summary <- mc_summary(data = data, mcnode = mcnode, 
                                 mc_name = node_name,
                                 keys_names = keys_names)
    }
    
    node_list[[i]][["summary"]] <- node_summary
    node_list[[i]][["keys"]] <- mc_summary_keys(node_summary)
  }
  
  if (!is.null(mcmodule)) {
    mcmodule$node_list <- node_list
    return(mcmodule)
  }
  
  return(node_list)
}