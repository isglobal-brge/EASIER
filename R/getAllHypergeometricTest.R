#' HyperGeometric Test for positions relative to Island
#'
#' HyperGeometric Test for positions relative to Island
#'
#' @param significative numerical. Variable to take in to account as significative variable, could be FDR, p-value,...
#' @param criteria dataframe with all relative positions to island values or chromatin status 15 or 18 Fetal Placenta
#' @param outputdir string. Output path to store file with results, by default results are written in current dir
#' @param outputfile string. File name to store results if no name is provided results are not written.
#'
#' @return Dataframe with hypergeometric test results for all criteria values
#'
#' @export
getAllHypergeometricTest <- function(significative, criteria, outputdir = ".", outputfile = NULL)
{
   positions <- unique(criteria)

   lregs <-  lapply(positions, function(x) getHypergeometricTest(significative, ifelse(criteria == x, "yes", "no"), x) )
   ans <- data.frame(matrix(unlist(lregs), nrow=length(positions), byrow=T))
   colnames(ans) <- names(lregs[[1]])
   rownames(ans) <- positions

   if(!is.null(outputfile)) {
      if(!is.null(outputdir) & !is.na(outputdir) & outputdir!='.')
         dir.create(outputdir, showWarnings = FALSE)

      # Output filename
      filename <- tools::file_path_sans_ext(basename(outputfile))

      write.csv(ans, paste0(file.path( outputdir),"/Hypergeomtest_",filename,".csv"))
   }

   return(ans)
}
