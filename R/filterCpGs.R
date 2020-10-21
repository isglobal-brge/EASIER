#' Filter all CpGs
#'
#' Filter all CpGs attending to multiple criteria values
#'
#' @param cpgs string vector with CpGs
#' @param tofilter vector with data to apply filter
#' @param filter vector with filter conditions
#'
#' @return string vector with CpGs not in filter criteria
#'
#' @export
filterCpGs <- function(cpgs, tofilter, filter)
{

   if(length(cpgs)!= length(tofilter))
      stop("Length differ between CpGs and criteria")

   fcrit <- grep(paste(filter,collapse = "|"), tofilter)

   filteredCpGs <- cpgs[-fcrit ]

   return(filteredCpGs)

}
