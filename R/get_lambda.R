#' Get lambda
#'
#' Get lambda
#'
#' @param cohort Dataframe with CpGs related data
#' @param cpval String with p-value column name
#'
#' @return lambda
#'
#' @export
get_lambda <- function( cohort, cpval )
{
   lambda <- qchisq(median(cohort[,cpval]), df = 1, lower.tail = FALSE) / qchisq(0.5, 1)

   return(lambda)
}
