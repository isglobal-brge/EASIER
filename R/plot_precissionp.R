#' Precission plot
#'
#' Precission plot
#'
#' @param x Data
#' @param filename optional, file name to write the plot with png format
#' @param main optional, string with plot title
#'
#' @return distribution plot
#'
#' @export
plot_precissionp <- function(x, filename, main= NULL, ...)
{

   if (is.null(main))
      main = "Precision Plot -  1/median(SE) vs sqrt(n)"

   p <- ggplot2::ggplot( data = x, mapping = ggplot2::aes( x = round(sqrt_N,2), y = round(invSE,2) ) ) +
      ggplot2::theme_bw() +
      ggplot2::geom_point( size = 3, ggplot2::aes( colour = cohort ) ) +
      #geom_line( aes( group = cohort , colour = cohort ) ) +
      ggplot2::ggtitle( main ) +
      ggplot2::theme( legend.position = "bottom",
                      legend.text = ggplot2::element_text(size=6),
                      legend.title = ggplot2::element_blank()) +
      ggplot2::labs( x = "sqrt(n)",
                     y = "inv SE") +
      ggplot2::scale_shape_manual( values = 0:7 )

   if(!is.null(filename))
      ggplot2::ggsave(filename,p)

   return(p)


}
