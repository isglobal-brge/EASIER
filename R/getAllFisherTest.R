#' Significative variable Regression for all relative positions to island
#'
#' Get regression by relative positions to island position taking in to account a significative variable parameter like FDR values. This function gets all chromatin state regressions and optionally store results in a file
#'
#' @param significative numerical. Variable to take in to account as significative variable, could be FDR, p-value,...
#' @param position dataframe with all relative positions to island values to perform regression
#' @param outputdir optional string. Output path to store file with results, by default results are written in current dir
#' @param outputfile optional string. File name to store results if no name is provided results are not written.
#' @param plots boolean. If plot is TRUE, plot results
#'
#' @return
#'
#' @export
getAllFisherTest <- function(significative, position, outputdir = ".", outputfile = NULL, plots = TRUE )
{

   positions <- getUniqueValues(position)

   lregs <-  lapply(positions, function(x) getFisherTest(significative, ifelse(position == x, "yes", "no"), x) )

   ans <- data.frame(matrix(unlist(lregs), nrow=length(positions), byrow=T))
   colnames(ans) <- names(lregs[[1]])

   if(!is.null(outputfile)) {
      if(!is.null(outputdir) & !is.na(outputdir) & outputdir!='.')
         dir.create(outputdir, showWarnings = FALSE)

      # Output filename
      filename <- tools::file_path_sans_ext(basename(outputfile))

      write.csv(ans, paste0(file.path( outputdir),"/Fisher_",filename,".csv"))

   }

   if( plots )
      plot_RelativetoIsland(ans, outputdir = outputdir ,outputfile = outputfile)


   return(ans)

}
