#' Depletion and Enrichment Plot
#'
#' Depletion and Enrichment Plot
#'
#' @param x Dataframe with Depletion and Enrichment data
#' @param outputdir string with relative path
#' @param main optional, string with title
#' @param xlab optional, string with xlab text
#'
#' @return distribution plot
#'
#' @export
plot_Depletion_Enrichment <- function(x, outputdir = '.', outputfile = NULL, main='', xlab='',...)
{

   if( !'Position' %in% colnames(x)) {
      x$Position <- rownames(x)
   }

   x$RelIsland <- factor(x[,"Position"], levels = (as.character(x[,"Position"])))
   x <- melt(x, id.vars = c("Position"), measure.vars = c("depletion", "enrichment"), value.name = 'OR' )


   nms <- names(x)
   x.plot <- nms['Position']

   p <- ggplot(x, aes(x = Position, y = OR, fill = factor(variable)  )) +
      geom_bar(stat="identity",  width = 0.5, position=position_dodge()) +
      theme_classic(base_size = 20) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_hline(yintercept = 1) +
      xlab(paste0("Relation to ",names(x)[1]))+
      ylab(paste0("prop.",names(x)[1]))

   if(!is.null(outputfile)) {
      if(!is.null(outputdir) & !is.na(outputdir) & outputdir!='.')
         dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)

      # Output filename
      filename <- tools::file_path_sans_ext(basename(outputfile))

      ggplot2::ggsave(paste0(file.path( outputdir),"/Decrept_Enrich_",filename,"_",names(x)[1],".pdf"), p)
   }

   return(p)
}
