#' Add cromatine information
#'
#' Add columns with crom15 information
#'
#' @param data dataframe used to add columns
#'
#' @return Original dataframe with cromatime related data
#'
#' @export
addCrom15Columns <- function(data)
{
   if(! "Name" %in% colnames(data))
      stop("CpG 'Name' field not found in data")

   # Adding phantom summary
   sub_crom15 <- crom15[,c(1,7:21)]

   data <- merge(data, sub_crom15, by.x = "Name", by.y = "HT12v4.ArrayAddress", sort = F)

   retur(data)

}
