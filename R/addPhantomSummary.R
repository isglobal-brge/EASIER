#' Add phantom summary
#'
#' Add summary from 'Phantom' annotated column in data
#'
#' @param data dataframe used to create columns
#'
#' @return Original dataframe with Phantom_S field data split in multiple column (TSS200, TSS1500, 5'UTR, 1stExon, Body, 3'UTR) with TRUE/FALSE
#'
#'
#' @export
addPhantomSummary <- function(data)
{
   if(! "Phantom" %in% colnames(data))
      stop("'Phantom' not found in data")

   # Adding phantom summary
   data$Phantom_S <- ifelse(grepl("low", data$Phantom), "low",
                            ifelse(grepl("high", data$Phantom),"high", ""))

   retur(data)

}
