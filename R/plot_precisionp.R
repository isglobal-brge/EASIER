#' Precision plot
#'
#' Precision plot
#'
#' @param x Data
#' @param filename optional, file name to write the plot with png format
#' @param main optional, string with plot title
#'
#' @return distribution plot
#'
#' @export
plot_precisionp <- function(x, filename, main= NULL, ...)
{

   if (is.null(main))
      main = "Precision Plot -  1/median(SE) vs sqrt(N)"

   p <- ggplot2::ggplot( data = x, mapping = ggplot2::aes( x = sqrt_N, y = invSE ) ) +
      ggplot2::theme_bw() +
      ggplot2::geom_point( size = 3, ggplot2::aes( colour = cohort ) ) +
      #geom_line( aes( group = cohort , colour = cohort ) ) +
      ggplot2::ggtitle( main ) +
      ggplot2::theme( legend.position = "bottom",
                      legend.text = ggplot2::element_text(size=6),
                      legend.title = ggplot2::element_blank()) +
      ggplot2::scale_shape_manual( values = 0:7 )

   if(!is.null(filename))
      ggplot2::ggsave(filename,p)

   return(p)


}
