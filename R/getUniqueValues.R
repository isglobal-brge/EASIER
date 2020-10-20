#' Get unique values in vector
#'
#' Get unique values from a Vector with strings separated by ";"
#'
#' @param x Vector with strings separated by ";"
#'
#' @return Vector with unique string values
#'
#' @export
getUniqueValues <- function(x)
{
   if(is.factor(x)){
      warning("x must be a character, conversion applied")
      x <- as.character(x)
   }

   temp <-  unlist(sapply(x, function(x) strsplit(x, ";")))
   ans <- unique(temp)
   ans <- ans[!is.na(ans)]

   return(ans)
}
