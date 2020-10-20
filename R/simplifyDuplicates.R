#' Get unique values in a cell
#'
#' Get unique values in a cell with data separated by ;
#'
#' @param x strings with duplicated data separated by ";"
#'
#' @return string with unique values separated by ;
#'
#' @export
simplifyDuplicates <- function(x)
{

   if(is.factor(x)){
      warning("x must be a character, conversion applied")
      x <- as.character(x)
   }

   tmp <- sapply(x, getUniqueValues)
   ans <- sapply(tmp, function(t) paste(t, collapse = ';'))

   return(ans)
}
