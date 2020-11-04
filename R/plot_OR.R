#' Plot OR
#'
#' Plot OR
#'
#' @param x Dataframe with OR data
#' @param outputdir string with relative path
#' @param main optional, string with title
#' @param xlab optional, string with xlab text
#'
#' @return distribution plot
#'
#' @export
plot_OR <- function(x, outputdir = '.', outputfile = NULL, main='', xlab='',...)
{

   x$RelIsland <- factor(x[,1], levels = (as.character(x[,1])))

   nms <- names(x)
   x.plot <- nms[1]

   # p <- ggplot(x, aes(x = OR, y = !!ensym(x.plot))) +
   #    geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") +
   #    geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height =
   #                      .2, color = "gray50") +
   #    geom_point(size = 3.5, color = "orange") +
   #    coord_trans(x = scales:::exp_trans(10)) +
   #    scale_x_continuous(breaks = log10(seq(0.1, 2.5, 0.1)), labels = seq(0.1, 2.5, 0.1),
   #                       limits = log10(c(0.09,2.5))) +
   #    theme_bw()+
   #    theme(panel.grid.minor = element_blank()) +
   #    ylab("") +
   #    xlab("Odds ratio") +
   #    annotate(geom = "text", y =1.1, x = log10(1.5),
   #             label = "Model p < 0.001\nPseudo R^2 = 0.10", size = 3.5, hjust = 0) +
   #    ggtitle("Feeding method and risk of obesity in cats")

   x[2:4] <- lapply(x[2:4], function(xf) as.numeric(as.character(xf)))

   ormin <- min(x$OR.inf)
   ormax <- max(x$OR.sup)

   p <- ggplot(x, aes(x = !!ensym(x.plot), y = OR)) +
      geom_bar(stat="identity", fill = "steelblue1", width = 0.5) +
      geom_errorbar(aes(ymin=OR.inf, ymax=OR.sup), width=0.2) +
      theme_classic(base_size = 20) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      xlab("") +
      ylab("Odds ratio")

   if(!is.null(outputfile)) {
      if(!is.null(outputdir) & !is.na(outputdir) & outputdir!='.')
         dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)

      # Output filename
      filename <- tools::file_path_sans_ext(basename(outputfile))

      ggplot2::ggsave(paste0(file.path( outputdir),"/OR_",filename,"_FDR_",names(x)[1],".pdf"), p)
   }

   return(p)
}
