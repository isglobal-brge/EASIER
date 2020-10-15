#' Chromosome State Plot
#'
#' Chromosome State plot
#'
#' @param x Dataframe with Chromosome State data
#' @param outputdir string with relative path
#' @param main optional, string with title
#' @param xlab optional, string with xlab text
#'
#' @return distribution plot
#'
#' @export
plot_chromosomestate <- function(x, outputdir = '.', outputfile = NULL, main='', xlab='',...)
{

   x$ChromStates <- factor(x$ChromStates, levels = (as.character(x$ChromStates)))
   p <- ggplot(x, aes(x = ChromStates, y = OR)) +
      geom_bar(stat="identity", fill = "steelblue1", width = 0.5) +
      geom_errorbar(aes(ymin=OR.inf, ymax=OR.sup), width=0.2) +
      theme_classic(base_size = 20) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_hline(yintercept = 1) +
      xlab("Chromatin States")

   if(!is.null(outputfile)) {
      if(!is.null(outputdir) & !is.na(outputdir) & outputdir!='.')
         dir.create(outputdir, showWarnings = FALSE)

      # Output filename
      filename <- tools::file_path_sans_ext(basename(outputfile))

      ggplot2::ggsave(paste0(file.path( outputdir),"/OR_",filename,".pdf"), p)
   }

   return(p)
}
