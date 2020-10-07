#' Write QC results to a file
#'
#' Write results from Quaity Control to a file. In filename the prefix _QCData is added in order to use this data in next steps like Meta-analysis
#'
#' @param cohort Dataframe with data
#' @param filename File name to write Quality Control results, func
#'
#' @return void
#'
#'
#' @export
write_QCData <- function(cohort, filename)
{
   # Write complete data to external file
   qc.fname.QCData <- paste0(filename, '_QCData.txt')
   suppressWarnings(
      write.table( as.data.frame(cohort[,c("probeID", "BETA", "SE", "P_VAL", "padj.fdr", "padj.bonf", "CpG_chrm", "CpG_beg", "CpG_end")]),
                   qc.fname.QCData, col.names = TRUE, row.names = FALSE, sep = '\t', append = FALSE, dec='.'))

}
