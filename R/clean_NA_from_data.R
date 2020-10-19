#' Remove NA from data
#'
#' Remove rows with NA data and return only complete cases.
#'
#' @param data dataframe. Data frame with na data
#' @return data frame without NAs
#'
#' @export
clean_NA_from_data <- function(data)
{
   return( data[complete.cases(data), ] )
}
