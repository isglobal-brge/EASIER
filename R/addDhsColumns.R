#' Add dhs information
#'
#' Add columns with dhs information
#'
#' @param data dataframe used to add dhs columns
#'
#' @return Original dataframe with dhs related data
#'
#'
#' @export
addDhsColumns <- function(data)
{
   if(! "Name" %in% colnames(data))
      stop("CpG 'Name' field not found in data")

   # Adding phantom summary
   sub_dhs <- dhs[,c(1,7:19)]

   data <- merge(data, sub_dhs, by.x = "Name", by.y = "HT12v4.ArrayAddress", sort = F)

   retur(data)

}
