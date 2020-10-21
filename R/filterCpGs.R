#' Filter all CpGs
#'
#' Filter all CpGs attending to multiple criteria values
#'
#' @param cpgs string vector with CpGs
#' @param position vector with data to apply filter
#' @param filter vector with filter conditions
#'
#' @return string vector with CpGs not in filter criteria
#'
#' @export
filterAllCpGs <- function(cpgs, criteria, filter)
{

   if(length(cpgs)!= length(criteria))
      stop("Length differ between CpGs and criteria")

   # Create dataframe with cpgs and criteria data
   data <- as.data.frame(cbind(cpgs, position))

   sapply(filter, filterCpGs, data = data )

   return(PMD.GRange)

}


#' Filter CpGs
#'
#' Filter all CpGs attending to a criteria values
#'
#' @param tofilter vector with data to filter
#' @param filter filter condition
#'
#' @return string vector with CpGs not in filter criteria
#'
#' @export
filterCpGs <- function(tofilter, filter)
{

   fdata <-  lapply(tofilter, function(tf) ifelse(tf == filter, "yes", "no") )

   return(fdata)

}
