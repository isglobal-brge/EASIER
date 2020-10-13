#' Print table with enrichment data
#'
#' Prints table with enrichment data
#'
#' @param x dataframe with data to print
#' @param columns columns to show in table
#'
#' @return void
#'
#'
#' @export
printTableEnrich <- function(x, columns=c(1:6)){
   kk <- as.data.frame(x)
   kk2 <- kk[kk$Count>2, columns]
   if (nrow(kk2)>=1)
      print(kable(kk2, row.names=FALSE))
   else
      cat("    No significant results \n")
   cat("\n")
}
