#' Call Rate filter
#'
#' filters out all CpGs with missing call rates exceeding the provided value
#'
#' @param data dataframe with data to be filtered
#' @param colname string, column name that contains number of samples for each CpG
#' @param pcmissing numeric, máximum percent of missing samples allowed
#' @param Ntotal numeric, Number of samples in cohort
#'
#' @return filtered dataframe if colname exists
#'
#' @export
filterLowRepresentedCpGsinCohort <- function(data, colname, pcmissing, Ntotal, fileresume)
{

   if( colname=='' | is.na(colname) | is.null(colname)) {
      return(data)
   }

   if( ! colname %in% names(data) ) {
      warning(paste0(colname, " does not exist in dataframe, no filter applied"))
      return(data)
   }

   totCpgs <- dim(data)[1]
   data$pecent <- data[colname]/Ntotal
   data <- data[which(data$pecent>=pcmissing),]
   remCpGs <- totCpgs - dim(data)[1]

   # Report descriptive exclussions to a descriptive file
   if(!is.null(fileresume)) {
      write(sprintf('\n# %s', strrep("-",36)), file = fileresume, append = TRUE)
      write(sprintf('# Remove CpGs with low representation: '), file = fileresume, append = TRUE)
      write(sprintf('# %s\n', strrep("-",36)), file = fileresume, append = TRUE)
      write(sprintf('# Minimum percentage : %s\n', pcmissing), file = fileresume, append = TRUE)
      write(sprintf('# Total CpGs in data : %d', dim(data)[1]), file = fileresume, append = TRUE)
      write(sprintf('# Number of excluded CpGs: %d', remCpGs), file = fileresume, append = TRUE)
      write(sprintf('# Total CpGs after exclusion : %d\n', dim(data)[1] ), file = fileresume, append = TRUE)
      write(sprintf('# Percent excluded CpGs: %f %%\n', ((remCpGs/totCpgs)*100 )), file = fileresume, append = TRUE)
   }

   return(data)

}

