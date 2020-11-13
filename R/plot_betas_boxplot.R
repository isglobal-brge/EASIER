#' Betas Boxplot
#'
#' Plot Betas boxplot
#'
#' @param x numeric list with data for each factor
#' @param filename, file name to write the boxplot with png format
#' @param main optional, string with plot title
#'
#' @return Betas boxplot
#'
#' @importFrom reshape2 melt
#' @import ggplot2 dplyr
#'
#' @export
plot_betas_boxplot <- function(x, filename= NULL, main= NULL, ...)
{

   if(is.null(filename))
      stop("No filename to store plot")

   if (is.null(main))
      main = "BETAS Boxplot"

   x <- melt(x)
   p <- ggplot2::ggplot( x, aes(x=factor(L1), y=value, color=factor(L1))) +
      ggplot2::ggtitle( main ) +
      ggplot2::geom_boxplot()+
      ggplot2::theme_bw() +
      ggplot2::theme( legend.position = "none",
                      axis.text.x = element_text( angle=90, hjust=1, vjust=0.5))

   ggplot2::ggsave(filename, p)

   return(p)

}
