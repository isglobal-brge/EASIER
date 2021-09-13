#' Gene Pisition Plots
#'
#' Gene Pisition plot
#'
#' @param x Dataframe with Gene Position data
#' @param outputdir string with relative path
#' @param main optional, string with title
#' @param xlab optional, string with xlab text
#'
#' @return distribution plot
#'
#' @export
plot_GenePosition <- function(x, outputdir = '.', outputfile = NULL, main='', xlab='',...)
{

   x$RelIsland <- factor(x[,1], levels = (as.character(x[,1])))
   # x[c("OR", "OR.inf","OR.sup","p-val")] <- lapply(x[c("OR", "OR.inf","OR.sup","p-val")], function(xf) as.numeric(levels(xf))[xf])
   x[c("OR", "OR.inf","OR.sup","p-val")] <- lapply(x[c("OR", "OR.inf","OR.sup","p-val")], function(xf) as.numeric(as.character(xf)))


   nms <- names(x)
   x.plot <- nms[1]

   p <- ggplot(x, aes(x = !!ensym(x.plot), y = OR)) +
      geom_bar(stat="identity", fill= "steelblue1", width = 0.5) +
      geom_errorbar(aes(ymin=OR.inf, ymax=OR.sup), width=0.2) +
      theme_classic(base_size = 20) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_hline(yintercept = 1) +
      xlab(paste0("Relation to ",names(x)[1]))

   if(!is.null(outputfile)) {
      if(!is.null(outputdir) & !is.na(outputdir) & outputdir!='.')
         dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)

      # Output filename
      filename <- tools::file_path_sans_ext(basename(outputfile))

      ggplot2::ggsave(paste0(file.path( outputdir),"/OR_",filename,"_",names(x)[1],".pdf"), p)
   }

   return(p)
}
