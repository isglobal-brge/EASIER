#' Significative variable Regression for all cromatin state
#'
#' Get regression by cromatin state taking in to account a significative variable parameter like FDR values. This function gets all chromatin state regressions and optionally store results in a file
#'
#' @param significative numerical. Variable to take in to account as significative variable, could be FDR, p-value,...
#' @param chromstates dataframe with all Cromatin state values to perform regression
#' @param outputdir string. Output path to store file with results, by default results are written in current dir
#' @param outputfile string. File name to store results if no name is provided results are not written.The suffix "RegressionFDR_States" is added to file name provided.
#' @param plots boolean. If plot is TRUE, plot results
#' @filename
#'
#' @return
#'
#' @export
getAllChromStateRegressions <- function(significative, chromstate, outputdir = ".", outputfile = NULL, plots = TRUE )
{

   lregs <-  lapply(colnames(chromstate), function(x) getChromStateRegression(significative, chromstate[,x], x) )

   ans <- data.frame(matrix(unlist(lregs), nrow=length(colnames(chromstate)), byrow=T))
   colnames(ans) <- names(lregs[[1]])

   if(!is.null(outputfile)) {
      if(!is.null(outputdir) & !is.na(outputdir) & outputdir!='.')
         dir.create(outputdir, showWarnings = FALSE)

      # Output filename
      filename <- tools::file_path_sans_ext(basename(outputfile))

      write.csv(ans, paste0(file.path( outputdir),"/OR_",filename,".csv"))

   }

   if( plots )
      plot_chromosomestate(ans, outputdir = outputdir ,outputfile = outputfile)


   return(ans)

}
