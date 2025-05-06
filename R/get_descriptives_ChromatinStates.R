#' Get descriptives Chromatine States
#'
#' Get counts, proportions and plots for Chromatine States CpGs attending to any criteria as Bonferroni, FDR
#'
#' @param data string. CpG
#' @param positions array string. with colnames containing chromatin state positions
#' @param criteria boolean. Boolean with criteria like Bonferroni, FDR,...
#' @param namecriteria string, criteria name for plot and to add to descriptives
#' @param outputdir string. Output path to store file with results, by default results are
#' written in current dir
#' @param outputfile string. File name to store results if no name is provided results are
#' not written.The suffix "_DescChromStates.txt" and "_PlotChromStates.pdf" are added
#' respectively for descriptive file and plots.
#'
#' @return File with descriptive and plots from CpG Annotation
#'
#'
#' @export
get_descriptives_ChromatineStates <- function(data, positions, criteria, namecriteria, outputdir = ".", outputfile = NULL)
{
   # Chromatine States for CpGs attending to criteria

   if( any( !positions %in% colnames(data)) ){
      stop( paste0( paste0(positions[which(!positins %in% data)], collapse = ","),
                    " column/s are not present in data"))
   }

   if( is.null(outputfile) ){
      stop("output file is needed to write descriptives")
   } else {
      if(!is.null(outputdir) & !is.na(outputdir) & outputdir!='.')
         dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
      # Output filename
      filename <- tools::file_path_sans_ext(basename(outputfile))
   }


   data <- data %>%
      mutate(ChromStates = pmap_chr(dplyr::select(., all_of(positions)), ~ {
         vals <- c(...)
         paste(positions[which(vals)], collapse = ", ")
      }))


   ctable <- as.data.frame.matrix(table(data$ChromStates, criteria))
   ptable <- as.data.frame.matrix(prop.table(ctable))

   desc.fname <- paste0(outputdir,"/",filename,"_",namecriteria,"_DescChromStates.txt")

   write(sprintf('\t\t============\n\t\t  Descriptive Resume for Chromatine States with %s criteria \n\t\t============\n',
                 namecriteria), file = desc.fname)
   write(sprintf('---------------------------------------------\n Model : %s\n---------------------------------------------\n',
                 filename), file = desc.fname, append = TRUE)
   write(sprintf('\nCount table : \n '), file = desc.fname, append = TRUE)
   suppressWarnings(write.table(ctable, file = desc.fname, append = TRUE))
   write(sprintf('\nProportion table : \n'), file = desc.fname, append = TRUE)
   suppressWarnings(write.table(ptable, file = desc.fname, append = TRUE))
   write(sprintf('\nChi-Square Test : \n'), file = desc.fname, append = TRUE)
   write(sprintf('\tX-squared : %f',chisq.test(ctable)[1]), file = desc.fname, append = TRUE)
   write(sprintf('\tdf : %f',chisq.test(ctable)[2]), file = desc.fname, append = TRUE)
   write(sprintf('\tp-value : %f',chisq.test(ctable)[3]), file = desc.fname, append = TRUE)

   pt <- as.data.frame(cbind(data$ChromStates, criteria))
   colnames(pt) <- c('ChromStates',namecriteria)
   ptc <- ggplot(pt, aes(ChromStates)) +
      geom_bar(aes(fill = criteria)) +
      scale_x_discrete(guide = guide_axis(angle = 90))

   plot.fname <- paste0(outputdir,"/",filename,"_",namecriteria,"_PlotChromStates.pdf")

   pdf(plot.fname)
      print(ptc)
   dev.off()

}
