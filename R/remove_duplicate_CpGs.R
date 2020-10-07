#' Remove duplicate CpGs
#'
#' Get a list with information related to duplicated values in data
#'
#' @param cohort Dataframe with data
#' @param descfile Optional, filename to write data description if NULL no descriptive values write
#' @param duplifile Optional, filename to write duplicates if NULL no results stored in a file
#'
#' @return cohort data without duplicated values
#'
#'
#' @export
remove_duplicate_CpGs <- function(cohort, colname, descfile = NULL, duplifile = NULL)
{

   # get CpG coulumn position in dataframe
   colIndex <- which(colnames(cohort)==colname)

   # Number of CpGs.
   ininCpG <- nrow(cohort)
   # if(!is.null(descfile)){
   #    write(sprintf('Descriptive information before remove duplicate CpGs \n'), file = descfile)
   #    write(sprintf('# Number of CpGs (Initially) : %d', ininCpG), file = descfile, append = TRUE)
   # }

   # Duplicates.
   dups <- which(duplicated(cohort[,colname]))
   #..# if(!is.null(descfile)) write(sprintf('# Number of duplicate CpGs: %d', length(dups)), file = qc.fname, append = TRUE)

   if (length(dups)) {
      warning(paste0('Removed ',length(dups) ,' duplicated probe names : ',paste(dups,', ')))

      # Report duplicates - only if there are duplicate CpGs
      if(!is.null(duplifile)) write.table(cohort[which(duplicated(cohort[,colname])),], duplifile, col.names = TRUE, row.names = TRUE, dec='.')

      cohort <- unique(cohort)

   }

   # Number of CpGs.
   nCpG <- nrow(cohort)

   # Summary - Previous inconsistent CpGs deletion
   summary <- apply(cohort[, -c( colIndex )], 2, summary)

   if(!is.null(descfile)){
      write(sprintf('Descriptive information before remove duplicate CpGs \n'), file = descfile)
      write(sprintf('# Number of CpGs (Initially) : %d', ininCpG), file = descfile, append = TRUE)
      write(sprintf('# Number of duplicate CpGs: %d', length(dups)), file = descfile, append = TRUE)
      write(sprintf('# Number of CpGs (after remove duplicates): %d \n', nCpG), file = descfile, append = TRUE)
      write(sprintf('# Data summary : \n'), file = descfile, append = TRUE)
      suppressWarnings( write.table(summary, file = descfile, append = TRUE, quote = TRUE, sep = '\t', dec='.') )
   }

   res <- cohort

}
