#' Match data or mc_nodes By Scenario or Groups
#' 
#' Match data tables or Monte Carlo nodes and with different dimensions by scenario or groups

#' Add Group IDs to Data Frames
#'
#' Creates a table summarizing origin farms with equivalent status
#'
#' @param data_x Primary dataframe
#' @param by Vector of column names to group by
#' @param data_y Optional secondary dataframe for matching
#' @return Dataframe or list of dataframes with added group IDs
#' @export
#' @examples
#' add_group_id(purchase_data, by=c("animal_category", "pathogen", "mov_id", "farm_id"))
add_group_id <- function(data_x, by, data_y=NULL) {
  if(!is.null(data_y)) {
    if(!all(by %in% names(data_x))) {
      stop(paste0(by[by %in% names(data_x)], " columns are not present in data_x"))
    }
    if(!all(by %in% names(data_y))) {
      stop(paste0(by[by %in% names(data_y)], " columns are not present in data_y"))
    }
    
    data_x$df <- "x"
    data_y$df <- "y"
    
    data_xy <- rbind(data_x[c(by,"df")], data_y[c(by,"df")]) %>%
      mutate(g_id=NULL,
             g_row=NULL) %>%
      group_by_at(vars(any_of(by))) %>%
      mutate(g_id=cur_group_id())
    
    data_x <- data_xy %>%
      filter(df=="x") %>%
      cbind(data_x[!names(data_x) %in% c(by,"df","g_id","g_row")]) %>%
      mutate(df=NULL,
             g_row=cur_group_rows()) %>%
      relocate(g_id, g_row) %>%
      ungroup()
    
    data_y <- data_xy %>%
      filter(df=="y") %>%
      cbind(data_y[!names(data_y) %in% c(by,"df","g_id","g_row")]) %>%
      mutate(df=NULL,
             g_row=cur_group_rows()) %>%
      relocate(g_id, g_row) %>%
      ungroup()
    
    return(list(x=data_x, y=data_y))
    
  } else {
    data_x <- data_x %>%
      group_by_at(vars(any_of(by))) %>%
      mutate(g_id=cur_group_id(),
             g_row=cur_group_rows()) %>%
      relocate(g_id,g_row) %>%
      ungroup()
    
    return(data_x)
  }
}

#' Match Groups Between Datasets
#'
#' @param x First dataset
#' @param y Second dataset
#' @param by Optional grouping variables
#' @return Matched dataset
#' @export
#' @examples
#' group_match(dataset1, dataset2, by=c("category", "group"))
group_match <- function(x, y, by=NULL) {
  if(is.null(by)) {
    message("Match by homogeneous groups (hg)")
    
    if(!"hg" %in% names(y) & !"hg" %in% names(x)) {
      stop("No homogeneous group index (hg) found and no match 'by' parameter specified")
    }
    
    if(!sum(x$scenario_id=="0")==sum(y$scenario_id=="0")) {
      stop("Differing number of rows in x and y. Specify variables to match in 'by' parameter.")
    }
    
    if(!"hg" %in% names(x)) {
      x$hg <- y$hg[y$scenario_id=="0"]
    }
    if(!"hg" %in% names(y)) {
      y$hg <- x$hg[x$scenario_id=="0"]
    }
  }
  
  y %>%
    add_group_id(by=by) %>%
    select(g_id, true_scenario=scenario_id, true_hg=hg) %>%
    left_join(add_group_id(x, by=by)) %>%
    mutate(scenario_id=true_scenario,
           hg=true_hg,
           true_hg=NULL,
           true_scenario=NULL) %>%
    relocate(hg, g_id, g_row, scenario_id)
}

#' Match Monte Carlo Nodes
#'
#' Matches monte carlo nodes with differing dimensions assuming the new variates 
#' are equal for scenario 0 for the smallest node
#'
#' @param mcmodule Monte Carlo module
#' @param mc_name Name of the Monte Carlo node
#' @param data Input data
#' @return Modified Monte Carlo node
#' @export
#' @examples
#' mc_match(mcmodule = purchase, mc_name="fattening_b_indir_contact_all", data=quarantine_data)
mc_match <- function(mcmodule, mc_name, data) {
  # Get previous node data
  prenode_hg <- mcmodule$node_list[[mc_name]][["hg"]]
  prenode_scenario <- mcmodule$node_list[[mc_name]][["scenario"]]
  match_prev_mcnode <- mcmodule$node_list[[mc_name]][["mcnode"]]
  
  if(is.null(match_prev_mcnode)) {
    stop(mc_name, "is null")
  }
  
  # Handle scalar to mcnode conversion
  if(!is.mcnode(match_prev_mcnode) & is.numeric(match_prev_mcnode)) {
    message(paste0(mc_name," to mcnode ", dim(match_prev_mcnode), "\n"))
    match_prev_mcnode <- mcdata(match_prev_mcnode, type="0", nvariates=length(prenode_hg))
  }
  
  prev_dim <- dim(match_prev_mcnode)
  
  # Get current data
  data_hg <- data$hg
  data_scenario <- data$scenario_id
  
  message("Matching dimensions by homogeneous groups provided (",
          max(prenode_hg,data_hg)," hg). From: ",
          length(unique(prenode_scenario)), " scenarios to ", 
          length(unique(data_scenario))," scenarios.")
  
  # Handle new scenarios
  new_scenario_hg <- data_hg[!data_scenario %in% prenode_scenario]
  
  for(i in new_scenario_hg) {
    mc_i <- extractvar(match_prev_mcnode,i)
    match_prev_mcnode <- addvar(match_prev_mcnode,mc_i)
  }
  
  new_dim <- dim(match_prev_mcnode)
  
  message(mc_name, " prev dim: [", paste(prev_dim, collapse=", "), 
          "], new dim: [", paste(new_dim, collapse=", "),"]")
  
  return(match_prev_mcnode)
}


#' Match Monte Carlo Nodes with Null Values
#'
#' Matches two mc_nodes with differing dimensions, handling null values
#'
#' @param mcmodule Monte Carlo module
#' @param mc_name_x First node name
#' @param mc_name_y Second node name
#' @param keys_names Names of key columns
#' @return List containing matched nodes and index
#' @export
#' @examples
#' mc_agg_match(mcmodule = intro, 
#'              mc_name_x="no_purchase_inf_agg", 
#'              mc_name_y="b_entry_agg")
mc_agg_match <- function(mcmodule, mc_name_x, mc_name_y, keys_names=c("pathogen")) {
  # Remove scenario_id from keys
  keys_names <- keys_names[!keys_names=="scenario_id"]
  
  # Get summaries and nodes
  summary_x <- mcmodule$node_list[[mc_name_x]][["summary"]]
  summary_y <- mcmodule$node_list[[mc_name_y]][["summary"]]
  
  summary_list <- add_group_id(summary_x, keys_names, summary_y)
  summary_x <- summary_list$x
  summary_y <- summary_list$y
  
  mcnode_x <- mcmodule$node_list[[mc_name_x]][["mcnode"]]
  mcnode_y <- mcmodule$node_list[[mc_name_y]][["mcnode"]]
  
  # Create cross join
  raw_cross_xy <- cross_join(
    summary_x[c("g_id","g_row", "scenario_id")],
    summary_y[c("g_id","g_row", "scenario_id")])
  
  # Find null groups
  g_null <- c(
    unique(raw_cross_xy$g_id.x[!raw_cross_xy$g_id.x %in% raw_cross_xy$g_id.y]),
    unique(raw_cross_xy$g_id.y[!raw_cross_xy$g_id.y %in% raw_cross_xy$g_id.x]))
  
  anti_filter_x <- raw_cross_xy$g_id.x %in% g_null
  anti_filter_y <- raw_cross_xy$g_id.y %in% g_null
  
  # Handle null scenarios
  anti_index_xy <- raw_cross_xy[anti_filter_x|anti_filter_y,] %>%
    filter(scenario_id.x==0|scenario_id.y==0) %>%
    filter((scenario_id.y==0 & !duplicated(paste(g_id.x, scenario_id.x))) |
             (scenario_id.x==0 & !duplicated(paste(g_id.y, scenario_id.y))) &
             !(scenario_id.y==0 & scenario_id.x==0 & 
                 (duplicated(paste(g_id.x, scenario_id.x))) |
                 duplicated(paste(g_id.y, scenario_id.y)))) %>%
    distinct()
  
  index_xy <- filter(raw_cross_xy, 
                     g_id.x==g_id.y & (scenario_id.x==0|scenario_id.y==0))
  
  index_xy <- bind_rows(index_xy,anti_index_xy) %>%
    distinct()
  
  # Match nodes
  null_x <- 0
  null_y <- 0
  
  # Process X node
  for(i in 1:nrow(index_xy)) {
    g_row_x_i <- index_xy$g_row.x[i]
    
    if(index_xy$g_id.y[i] %in% summary_x$g_id) {
      mc_i <- extractvar(mcnode_x,g_row_x_i)
    } else {
      mc_i <- extractvar(mcnode_x, 1)-extractvar(mcnode_x, 1)
      null_x <- null_x+1
    }
    
    if(!exists("mcnode_x_match")) {
      mcnode_x_match <- mc_i
    } else {
      mcnode_x_match <- addvar(mcnode_x_match,mc_i)
    }
  }
  
  # Process Y node
  for(i in 1:nrow(index_xy)) {
    g_row_y_i <- index_xy$g_row.y[i]
    
    if(index_xy$g_id.x[i] %in% summary_y$g_id) {
      mc_i <- extractvar(mcnode_y,g_row_y_i)
    } else {
      mc_i <- extractvar(mcnode_y, 1)-extractvar(mcnode_y, 1)
      null_y <- null_y+1
    }
    
    if(!exists("mcnode_y_match")) {
      mcnode_y_match <- mc_i
    } else {
      mcnode_y_match <- addvar(mcnode_y_match,mc_i)
    }
  }
  
  # Log results
  message(mc_name_x, " prev dim: [", paste(dim(mcnode_x), collapse=", "),
          "], new dim: [", paste(dim(mcnode_x_match), collapse=", "),
          "], ", null_x, " null matches")
  
  message(mc_name_y, " prev dim: [", paste(dim(mcnode_y), collapse=", "),
          "], new dim: [", paste(dim(mcnode_y_match), collapse=", "),
          "], ", null_y, " null matches")
  
  # Create index
  index <- summary_x[ifelse(index_xy$g_id.x %in% summary_x$g_id,
                            index_xy$g_row.x,
                            index_xy$g_row.y),
                     c(keys_names)]
  
  scenario <- ifelse(
    summary_x[index_xy$g_row.x,]$scenario_id=="0",
    ifelse(summary_y[index_xy$g_row.y,]$scenario_id=="0","0",
           summary_y[index_xy$g_row.y,]$scenario_id),
    summary_x[index_xy$g_row.x,]$scenario_id)
  
  index$scenario_id <- scenario
  
  # Return results
  result <- list(mcnode_x_match, mcnode_y_match, index)
  names(result) <- c(paste0(mc_name_x,"_match"),
                     paste0(mc_name_y,"_match"),
                     "index")
  
  return(result)
}