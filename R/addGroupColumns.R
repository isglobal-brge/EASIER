#' Add group columns from annotation data
#'
#' Add group columns from UCSC_RefGene_Group annotated column. UCSC_RefGene_Group column must exists on inpu dataframe
#'
#' @param data dataframe used to create columns
#'
#' @return Original dataframe with UCSC_RefGene_Group field data split in multiple column (TSS200, TSS1500, 5'UTR, 1stExon, Body, 3'UTR) with TRUE/FALSE
#'
#'
#' @export
addGroupColumns <- function(data)
{
   if(! "UCSC_RefGene_Group" %in% colnames(data))
      stop("'UCSC_RefGene_Group' not found in data")

   # Adding group columns
   data$TSS200 <- ifelse(grepl("TSS200", data$UCSC_RefGene_Group), T, F)
   data$TSS1500 <- ifelse(grepl("TSS1500", data$UCSC_RefGene_Group), T, F)
   data$UTR5 <- ifelse(grepl("5'UTR", data$UCSC_RefGene_Group), T, F)
   data$FirstExon <- ifelse(grepl("1stExon", data$UCSC_RefGene_Group), T, F)
   data$Body <- ifelse(grepl("Body", data$UCSC_RefGene_Group), T, F)
   data$UTR3 <- ifelse(grepl("3'UTR", data$UCSC_RefGene_Group), T, F)

   retur(data)

}
