#' Distribution plot
#'
#' Distribution plot
#'
#' @param x Data
#' @param main optional, string with title
#' @param xlab optional, string with xlab text
#'
#' @return distribution plot
#'
#' @export
plot_distribution <- function(x, main='', xlab='',...)
{

   h <- hist(x, breaks = 100, plot = FALSE)
   d <- density(x)
   plot(h, freq = FALSE, col = gg_colors(2)[2], border = 'white',
        main = main, xlab = xlab, ylim = c(0, max(d$y, h$density)))
   lines(d, col = gg_colors(2)[1], lwd = 2)

}
