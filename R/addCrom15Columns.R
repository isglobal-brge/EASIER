#' Add Chromatine states information
#'
#' Add columns with Chromatine states information
#'
#' @param data dataframe used to add columns
#' @param column numeric index or column name column with CpG id.
#'
#' @return Original dataframe with cromatime related data
#'
#' @export
addCrom15Columns <- function(data, cpgcol)
{
   data("crom15")
   # Adding phantom summary
   sub_crom15 <- crom15[,c(2,8:22)]

   data <- merge(data, sub_crom15, by.x = cpgcol, by.y = "HT12v4.ArrayAddress", sort = F)

   return(data)

}
