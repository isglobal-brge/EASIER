#' Get descriptives gene position
#'
#' Get counts, proportions and plots for gene position attending to any criteria as Bonferroni, FDR
#'
#' @param position string. CpG Gene position (Body, TSS1500, TSS200, 3'UTR, 1stExon, 5'UTR )
#' @param criteria boolean. Boolean with criteria like Bonferroni, FDR,...
#' @param namecriteria string, criteria name for plot and to add to descriptives
#' @param outputdir string. Output path to store file with results, by default results are written in current dir
#' @param outputfile string. File name to store results if no name is provided results are not written.The suffix "_DescGene.txt" and "_PlotGene.pdf" are added respectively for descriptive file and plots.
#'
#' @return File with descriptive and plots from CpG Annotation
#'
#'
#' @export
get_descriptives_GenePosition <- function(position, criteria, namecriteria, outputdir = ".", outputfile = NULL)
{
   # Gene Positions for CpGs attending to criteria

   if( is.null(outputfile) ){
      stop("output file is needed to write descriptives")
   }else {
      if(!is.null(outputdir) & !is.na(outputdir) & outputdir!='.')
         dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
      # Output filename
      filename <- tools::file_path_sans_ext(basename(outputfile))
   }

   sposition <- simplifyDuplicates(position)


   ctable <- as.data.frame.matrix(table(sposition, criteria))
   if(dim(ctable)[2]>0) {
      ptable <- as.data.frame.matrix(prop.table(ctable))
   }

   desc.fname <- paste0(outputdir,"/",filename,"_",namecriteria,"_DescGene.txt")

   write(sprintf('\t\t\t\t============\n\t\t\t\t  Descriptive Resume for Gene Position with %s criteria \n\t\t\t\t============\n', namecriteria), file = desc.fname)
   write(sprintf('-------------------\n Model : %s\n-------------------\n',filename), file = desc.fname, append = TRUE)
   write(sprintf('\nCount table : \n '), file = desc.fname, append = TRUE)
   suppressWarnings(write.table(ctable, file = desc.fname, append = TRUE))

   if( exists("ptable")){
      write(sprintf('\nProportion table : \n'), file = desc.fname, append = TRUE)
      suppressWarnings(write.table(ptable, file = desc.fname, append = TRUE))
   }

   if(dim(ctable)[2]>0)
   {
      write(sprintf('\nChi-Square Test : \n'), file = desc.fname, append = TRUE)
      write(sprintf('\tX-squared : %f',chisq.test(ctable)[1]), file = desc.fname, append = TRUE)
      write(sprintf('\tdf : %f',chisq.test(ctable)[2]), file = desc.fname, append = TRUE)
      write(sprintf('\tp-value : %f',chisq.test(ctable)[3]), file = desc.fname, append = TRUE)
      pt <- as.data.frame(cbind(sposition, criteria))
      colnames(pt) <- c('position', namecriteria)
      ptc <- ggplot(pt, aes(sposition)) +
         geom_bar(aes(fill = criteria)) +
         scale_x_discrete(guide = guide_axis(angle = 90))

      plot.fname <- paste0(outputdir,"/",filename,"_",namecriteria,"_PlotGene.pdf")
      ggplot2::ggsave(plot.fname, ptc, device = pdf)
   } else {
      write(sprintf('\nChi-Square Test : \n'), file = desc.fname, append = TRUE)
      write(sprintf('\tNo test to perform '), file = desc.fname, append = TRUE)
   }


}
