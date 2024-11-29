#' Adapted Version of Function 'signif'
#' 
#' This function adapts base-function \code{\link{signif}}
#' by always returning integer values in case the number of
#' requested significant digits is less than the the number of
#' digits in front of the decimal separator.
#' 
#' @param x			(numeric) value to be rounded to the desired number
#' 					of significant digits
#' @param digits	(integer) number of significant digits
#' @param force		(logical) TRUE = force the return value to have at least 4 significant
#' 					digits, i.e. to integers with less digits zeros will be appended after
#' 					the decimal separator, otherwise the return value will be casted from
#' 					character to numeric
#' @param ...		additional parameters
#' 
#' @return 	number with 'digits' significant digits, if 'force=TRUE' "character" objects will be
#' 			returned otherwise objects of mode "numeric"
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}


Signif <- function(x, digits=4, force=TRUE, ...)
{
  call 	<- match.call()
  manyX 	<- call$manyX
  if(is.null(manyX))
    manyX <- FALSE
  stopifnot(is.numeric(x))
  if(length(x) > 1)
    return(sapply(x, Signif, digits=digits, manyX=TRUE))
  
  if(!manyX && "coef.gnm" %in% class(x))		# assign name to single gnm-coefficient
  {
    x <- as.numeric(x)
    names(x) <- "beta1"
  }
  Ndbc	<- nchar(substr(as.character(x), 1, regexpr("\\.", as.character(x))-1))
  x 		<- signif(x, ifelse(Ndbc > digits, Ndbc, digits))
  NcX		<- nchar(x)
  comma	<- grepl("\\.", x)
  if(comma)
    NcX <- NcX - 1
  if(NcX < digits)
    x <- paste0(x, ifelse(comma, "", "."), paste(rep(0, digits-NcX), collapse=""))
  if(!force)
    x <- as.numeric(x)
  x
}
