#' Test duplicate CpGs
#'
#' Get a list with information related to duplicated values in data
#'
#' @param cohort Dataframe with data
#' @param colname string, column name to search duplicates
#' @param descfile Optional, filename to write data description if NULL no descriptive values write
#'
#' @return cohort data without duplicated values
#'
#' @export
test_duplicate_CpGs <- function(cohort, colname,  duplifile = NULL)
{

   # get CpG coulumn position in dataframe
   colIndex <- which(colnames(cohort)==colname)

   # Number of CpGs.
   ininCpG <- nrow(cohort)

   # Duplicates.
   dups <- which(duplicated(cohort[,colname]))
   #..# if(!is.null(descfile)) write(sprintf('# Number of duplicate CpGs: %d', length(dups)), file = qc.fname, append = TRUE)

   if (length(dups)>0) {
      # Report duplicates - only if there are duplicate CpGs
      if(!is.null(duplifile)) write.table(cohort[which(duplicated(cohort[,colname])),], duplifile, col.names = TRUE, row.names = TRUE, dec='.')

      stop(paste0(length(dups) ,' duplicates found with probe names : ',paste(cohort[which(duplicated(cohort[,colname])),colname],', ')))
   }


}
