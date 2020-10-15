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

   if( is.null(outputfile) ){
      stop("output file is needed to write descriptives")
   }else {
      if(!is.null(outputdir) & !is.na(outputdir) & outputdir!='.')
         dir.create(outputdir, showWarnings = FALSE)
      # Output filename
      filename <- tools::file_path_sans_ext(basename(outputfile))
   }

   ctable <- as.data.frame.matrix(table(position, criteria))
   ptable <- as.data.frame.matrix(prop.table(ctable))

   desc.fname <- paste0(outputdir,"/",filename,"_",namecriteria,"_DescIslands.txt")

   write(sprintf('\t\t\t\t============\n\t\t\t\t  Descriptive Resume for Relative position to Island with %s criteria \n\t\t\t\t============\n', namecriteria), file = desc.fname)
   write(sprintf('-------------------\n Model : %s\n-------------------\n',filename), file = desc.fname, append = TRUE)
   write(sprintf('\nCount table : \n '), file = desc.fname, append = TRUE)
   suppressWarnings(write.table(ctable, file = desc.fname, append = TRUE))
   write(sprintf('\nProportion table : \n'), file = desc.fname, append = TRUE)
   suppressWarnings(write.table(ptable, file = desc.fname, append = TRUE))
   write(sprintf('\nChi-Square Test : \n'), file = desc.fname, append = TRUE)
   write(sprintf('\tX-squared : %f',chisq.test(ctable)[1]), file = desc.fname, append = TRUE)
   write(sprintf('\tdf : %f',chisq.test(ctable)[2]), file = desc.fname, append = TRUE)
   write(sprintf('\tp-value : %f',chisq.test(ctable)[3]), file = desc.fname, append = TRUE)

   pt <- as.data.frame(cbind(position, criteria))
   colnames(pt) <- c('position',namecriteria)
   ptc <- ggplot(pt, aes(position)) +
      geom_bar(aes(fill = criteria))

   plot.fname <- paste0(outputdir,"/",filename,"_",namecriteria,"_PlotIslands.pdf")
   ggplot2::ggsave(plot.fname, ptc, device = pdf)

}
