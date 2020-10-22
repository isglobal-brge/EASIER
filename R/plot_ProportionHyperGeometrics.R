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
plot_ProportionHyperGeometrics <- function(x, outputdir = '.', outputfile = NULL, main='', xlab='',...)
{

   p<-ggplot(x,aes(x = Criteria,y = value)) +
      geom_bar(aes(fill = variable),stat = "identity",position = "dodge") +
      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
      scale_fill_manual("Legend", values = c("All_CpGs" = "grey45", "Sig_CpGs" = "green3", "Sig_Hypo_CpGs" = "red3", "Sig_Hyper_CpGs" = "blue3"))  +
      labs(y="Proportion")

   print(p)
   if(!is.null(outputfile)) {
      if(!is.null(outputdir) & !is.na(outputdir) & outputdir!='.')
         dir.create(outputdir, showWarnings = FALSE)

      # Output filename
      filename <- paste0(file.path( outputdir),"/Prop_",tools::file_path_sans_ext(basename(outputfile)),"_States_FP.pdf")

      ggsave(filename)
         p
      dev.off()

   }

   print(p)

   return(p)

}
