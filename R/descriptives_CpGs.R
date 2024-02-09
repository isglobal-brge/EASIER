#' CpGs descriptive data
#'
#' Get CpG descriptive data and optionally, writes a file with descriptive data from CpGs
#'
#' @param cohort Dataframe with CpGs data
#' @param columns Integer vector with column position that we are interested in
#' @param filename Optional, filename to write results, if NULL no results are writed in a file
#' @param nsamples numeric, Optional, samples in cohort
#' @param before bolean, Optional, if true, prints text "Before QC" , else, prints "After QC"
#'
#' @return list with descriptive values
#' \itemize{
#'   \item{nCpGs}{Number of CpGs in data}
#'   \item{BETA}{Descriptive data for beta value}
#'   \item{SE}{Descriptive data for  standard Error}
#'   \item{P_VAL}{Descriptive data for p-value}
#' }
#'
#'
#' @export
descriptives_CpGs <- function(cohort, columns, filename = NULL, nsamples = NULL, artype = '450K', before = TRUE)
{
   #
   desc <- c('min.', '1Q', 'Median', 'Mean', '3Q.', 'max')

   # Number of CpGs.
   nCpG <- nrow(cohort)

   # Summary - Previous inconsistent CpGs deletion
   summary <- apply(cohort[, columns], 2, summary)

   if(!is.null(filename)){
      write(sprintf('\n# %s', strrep("-",29)), file = filename, append = TRUE)
      if(before == TRUE){
         write(sprintf('# Summary before QC : '), file = filename, append = TRUE)
      }else{
         write(sprintf('Summary after QC : '), file = filename, append = TRUE)
      }
      write(sprintf('# %s \n', strrep("-",29)), file = filename, append = TRUE)

      if(before == TRUE){
         # merge with array type (useful for EPIC arrays used as 450K arrays in analysis)
         write(sprintf('# Number of CpGs (originally): %d \n', nCpG), file = filename, append = TRUE)
         commonCpGs <- restrict_CpGs_to_artype(cohort$probeID, artype)
         cohort <- cohort[which(cohort$probeID %in% commonCpGs),]
         nCpG <- nrow(cohort)
         write(sprintf('# Number of CpGs (after match with array type): %d \n', nCpG), file = filename, append = TRUE)
      } else {
         # Resume File (with all QC analyses)
         summaryResFileName <- str_flatten(str_split_1(filename, "/")[1:(length(str_split_1(filename, "/") )-2)], "/")
         summaryResFileName <- paste(summaryResFileName, "tmp_postQC.txt", sep="/")

         # Write descriptives
         write(sprintf('# Number of CpGs: %d \n', nCpG), file = filename, append = TRUE)

         summAll <- as.data.frame(cbind( artype, nsamples, nCpG,  t(summary[,1]), t(summary[,2]), t(summary[,3]) ))
         colnames(summAll) <- c("array", "NSamples", "N Cpgs after QC", paste0("Coef_",desc), paste0("SE_",desc), paste0("pval_",desc) )

         # Write full summary file
         if(!file.exists(summaryResFileName)) {
            write.table(summAll, summaryResFileName,col.names = TRUE, append = FALSE, row.names = FALSE, sep="\t")
         }else {
            write.table(summAll, summaryResFileName, col.names = FALSE, append = TRUE, row.names = FALSE, sep="\t")
         }
      }

      if(!is.null(nsamples)){
         write(sprintf('# Number of samples: %d \n', nsamples), file = filename, append = TRUE)
      }
      write(sprintf('# Descriptive : \n'), file = filename, append = TRUE)
      suppressWarnings( write.table(summary, file = filename, append = TRUE, quote = TRUE, sep = '\t', dec='.') )
   }

   if( before == TRUE) {
      return(cohort)
   }

}


