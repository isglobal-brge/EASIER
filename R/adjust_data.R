#' Adjust p-value by Bonferroni or FDR
#'
#' Update data with Bonferroni and FDR adjustment
#'
#' @param cohort Dataframe with CpGs related data
#' @param cpval String with p-value column name
#' @param bn Boolean if TRUE appends bonferroni adjustment column
#' @param fdr Boolean if TRUE appends FDR adjustment column
#' @param filename Optional, filename to write results, if NULL no results are writed in a file
#' @param N Optional, informative value to write to a file with other descriptive data related to significative p-values.
#'
#' @return Data frame with original data with bonferroni and fdr adjusted data sorted by p-value
#'
#' @export
adjust_data <- function( cohort, cpval, bn = TRUE, fdr = TRUE, filename = NULL, N = 0)
{

   # Resume File (with all QC analyses)
   summaryResFileName <- str_flatten(str_split_1(filename, "/")[1:(length(str_split_1(filename, "/") )-2)], "/")
   summaryResFileName <- paste(summaryResFileName, "tmp_postQCAdj.txt", sep="/")

   cohort$padj.bonf <-  ifelse( cohort[,cpval]<0.05/dim(cohort)[1], 'yes','no')
   cohort$padj.fdr <- p.adjust(cohort[,cpval], method = "fdr")
   lambda <- qchisq(median(cohort[,cpval]), df = 1, lower.tail = FALSE) / qchisq(0.5, 1)

   # Sort data by p-val
   cohort <- cohort[order(cohort[,cpval]),]

   # Write resume to external file
   if(!is.null(filename)) {
      write(sprintf('\n# %s', strrep("-",21)), file = filename, append = TRUE)
      write(sprintf('# Significative CpGs : '), file = filename, append = TRUE)
      write(sprintf('# %s \n', strrep("-",21)), file = filename, append = TRUE)

      write(sprintf('# With significative p-value : %s',sum(cohort[,cpval] < 0.05) ), file = filename, append = TRUE)
      write(sprintf('# With significative p-value adjusted by FDR : %s', sum(cohort$padj.fdr<0.05)), file = filename, append = TRUE)
      write(sprintf('# With significative p-value adjusted by Bonferroni : %s',sum(cohort$padj.bonf =='yes') ), file = filename, append = TRUE)
      write(sprintf('# N : %s', N ), file = filename, append = TRUE)
      write(sprintf('# lambda : %s', lambda ), file = filename, append = TRUE)
   }

   summAll <- as.data.frame(cbind(lambda, sum(cohort[,cpval] < 0.05) ,  sum(cohort$padj.fdr<0.05), sum(cohort$padj.bonf =='yes')))
   names(summAll) <- c('lambda', 'N_Sig_Nominal', 'N_Sig_FDR', 'N_Sig_BN')

   # Write full summary file
   if(!file.exists(summaryResFileName)) {
      write.table(summAll, summaryResFileName,col.names = TRUE, append = FALSE, row.names = FALSE, sep = "\t" )
   }else {
      write.table(summAll, summaryResFileName, col.names = FALSE, append = TRUE, row.names = FALSE, sep = "\t" )
   }


   return(cohort)
}
