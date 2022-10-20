#' Betas Boxplot
#'
#' Plot Betas boxplot
#'
#' @param x numeric list with data for each factor
#' @param filename, file name to write the boxplot with png format
#' @param main optional, string with plot title
#' @param outliers optional, boolean by default = FALSE, if true,
#' plot_betas_boxplot function plots outliers
#'
#' @return Betas boxplot
#'
#' @importFrom reshape2 melt
#' @import ggplot2 dplyr
#'
#' @export
plot_betas_boxplot <- function(x, filename= NULL, main= NULL, outliers = FALSE, ...)
{

   if(is.null(filename))
      stop("No filename to store plot")

   if (is.null(main))
      main = "BETAS Boxplot"

   x <- melt(x)

   if(outliers == FALSE){
      p <- ggplot2::ggplot( x, aes(x=factor(L1), y=value, color=factor(L1))) +
         ggplot2::ggtitle( main ) +
         ggplot2::geom_boxplot(outlier.shape = NA) +
         ggplot2::theme_bw() +
         ggplot2::theme( legend.position = "none",
                         axis.text.x = element_text( angle=90, hjust=1, vjust=0.5))

      # compute lower and upper whiskers
      ylim = boxplot.stats(x$value)$stats[c(1, 5)]

      # scale y limits based on ylim1
      pf = p + coord_cartesian(ylim = ylim*1.05)

   } else {
      pf <- ggplot2::ggplot( x, aes(x=factor(L1), y=value, color=factor(L1))) +
         ggplot2::ggtitle( main ) +
         ggplot2::geom_boxplot() +
         ggplot2::theme_bw() +
         ggplot2::theme( legend.position = "none",
                         axis.text.x = element_text( angle=90, hjust=1, vjust=0.5))
   }

   png( filename, type = "cairo")
      print(pf)
   dev.off()
   # ggplot2::ggsave(filename, pf)

   return(pf)

}
