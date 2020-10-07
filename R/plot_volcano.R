#' Volcano plot
#'
#' Volcano plot
#'
#' @param x Data
#' @param cbeta string with beta column name
#' @param cpval string with p-value column name
#' @param main optional, string with title
#' @param xlab optional, string with xlab text
#'
#' @return distribution plot
#'
#' @export
plot_volcano <- function(x, cbeta, cpval, nCpGs, main='',...)
{

   bt <- 0.01
   pt <- 3
   nCpGs <- nrow(x)

   colors <- rep(gray(0.75, 0.25), nCpGs)
   colors[ x[,cbeta] >  bt & -log10(x[,cpval]) > 3] <- gg_colors(2, 0.5)[1] # Red
   colors[x[,cbeta] < -bt & -log10(x[,cpval]) > 3] <- gg_colors(2, 0.5)[2] # Blue
   plot(x[,cbeta], -log10(x[,cpval]), col = colors,
        main = main, xlab = 'Beta', ylab = '-log10 p-value')
   abline(h = pt, v = c(-bt, bt), lty = 'dotted')

}


# Custom graphical options.
gg_colors <- function(n, alpha = NULL) {
   hues = seq(15, 375, length = n + 1)
   hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}
