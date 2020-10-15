#' Get descriptives relatives to Island
#'
#' Get counts, proportions and plots for relative position to Island CpGs attending to any criteria as Bonferroni, FDR
#'
#' @param position string. CpG relative position to Island (Island, N_Shelf, N_Shore, OpenSea, S_Shelf, S_Shore )
#' @param criteria boolean. Boolean with criteria like Bonferroni, FDR,...
#' @param namecriteria string, criteria name for plot and to add to descriptives
#' @param outputdir string. Output path to store file with results, by default results are written in current dir
#' @param outputfile string. File name to store results if no name is provided results are not written.The suffix "_DescIslands.txt" and "_PlotIslands.pdf" are added respectively for descriptive file and plots.
#'
#' @return File with descriptive and plots from CpG Annotation
#'
#'
#' @export
get_descriptives_RelativetoIsland <- function(position, criteria, namecriteria, outputdir = ".", outputfile = NULL)
{
   # Relation to Island for CpGs attending to criteria

   if( is.null(outputfile) )
      stop("output file is needed to write descriptives")
   else {
      if(!is.null(outputdir) & !is.na(outputdir) & outputdir!='.')
         dir.create(outputdir, showWarnings = FALSE)
      # Output filename
      filename <- tools::file_path_sans_ext(basename(outputfile))
   }


   ctable <- as.data.frame.matrix(table(position, criteria))
   ptable <- as.data.frame.matrix(prop.table(position, criteria))

   desc.fname <- paste0(filename,"_",criteria,"_DescIslands.txt")
   print(paste0("Output file : ",qc.fname))

   write(sprintf('\t\t\t\t============\n\t\t\t\t  Descriptive Resume for Relative position to Island with %s criteria \n\t\t\t\t============\n', criteria), file = desc.fname)
   write(sprintf('-------------------\n Model : %s\n-------------------\n',metaname), file = desc.fname, append = TRUE)
   write(sprintf('\nCount table : \n '), file = desc.fname, append = TRUE)
   write(ctable, file = desc.fname, append = TRUE)
   print("1")
   write(sprintf('\nProportion table : \n'), file = desc.fname, append = TRUE)
   write(ptable, file = desc.fname, append = TRUE)
   print("2")
   write(sprintf('\nChi-Square Test : \n'), file = desc.fname, append = TRUE)
   write(chisq.test(ctable), file = desc.fname, append = TRUE)
   print("3")


   pt <- as.data.frame(cbind(position, criteria))
   colnames(pt) <- c('position',namecriteria)
   ptc <- ggplot(pt, aes(position)) + geom_bar(aes(fill = criteria))

   plot.fname <- paste0(filename,"_",criteria,"_PlotIslands.pdf")
   ggplot2::ggsave(plot.fname, ptc)



}
