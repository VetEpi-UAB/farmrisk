#' Replace NA and Infinite Values in mcnode Objects
#'
#' This function replaces NA and infinite values in mcnode objects with a specified value.
#'
#' @param mcnode An mcnode object containing NA or infinite values
#' @param na_value Numeric value to replace NA and infinite values (default = 0)
#'
#' @return An mcnode object with NA and infinite values replaced by na_value
#'
#' @examples
#' # Create a sample mcnode with NA values
#' sample_mcnode <- mcnode(matrix(c(1, NA, 3, Inf, 5), nrow = 5))
#' # Replace NA and Inf with 0
#' clean_mcnode <- mcnode_na_rm(sample_mcnode)
#'
#' @export
#'
#' @seealso \code{\link[mc2d]{is.na.mcnode}}
mcnode_na_rm <- function(mcnode, na_value = 0) {
  replace(mcnode, is.na(mcnode) | is.infinite(mcnode), na_value)
}