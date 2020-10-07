#' CpGs descriptive data
#'
#' Get CpG descriptive data and optionally, writes a file with descriptive data from CpGs
#'
#' @param cohort Dataframe with CpGs data
#' @param columns Integer vector with column position that we are interested in
#' @param filename Optional, filename to write results, if NULL no results are writed in a file
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
descriptives_CpGs <- function(cohort, columns, filename = NULL)
{
   # Number of CpGs.
   nCpG <- nrow(cohort)

   # Summary - Previous inconsistent CpGs deletion
   summary <- apply(cohort[, columns], 2, summary)

   if(!is.null(filename)){
      write(sprintf('Summary for CpGs \n'), file = filename)
      write(sprintf('# Number of CpGs: %d \n', nCpG), file = filename, append = TRUE)
      write(sprintf('# Data summary : \n'), file = filename, append = TRUE)
      suppressWarnings( write.table(summary, file = filename, append = TRUE, quote = TRUE, sep = '\t', dec='.') )
   }

   res <- list("nCpGs" = nCpG,
               "descriptives" = summary )

}


