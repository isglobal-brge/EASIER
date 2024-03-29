#' Plot enrichment summary
#'
#' Plot enrichment summary from GO and KEGG
#'
#' @param df Dataframe with enrichment data
#' @param xl string with x axis title
#' @param yl string with y axis title
#' @param tit string with main title
#' @param filename optional string. File name to store results if no name is provided results are not written.
#'
#' @return plot
#'
#' @export
plot_Summary <- function(df, xl, yl, tit, filename = NULL)
{
   df <- df[df$N>10,]
   mask <- sapply(df$TERM, nchar)>50
   df$TERM[mask] <- paste0(substr(df$TERM[mask], 1, 45),
                                  "_trunc")
   p <- ggplot(df, aes(x = ONTOLOGY, y = TERM, col = FDR)) +
      geom_point(aes(size = N)) +
      xlab(xl) + ylab(yl) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ggtitle(tit) + scale_size_continuous(limits=c(1,40)) +
      scale_colour_gradient(trans="reverse", low = "light green",
                            high = "dark green", na.value="transparent",
                            limits=c(0.05,0.0001))

   if(!is.null(filename)){
      #..# ggplot2::ggsave(filename, p)
      png( filename, type = "cairo")
         print(p)
      dev.off()
   }

}


