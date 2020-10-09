#' Get unique genes
#'
#' Get unique genes from annotated results
#'
#' @param x Annotated genes vector
#' @param entrez boolean. Optional if TRUE (default option), use org.Hs.egSYMBOL2EG to return entrez id
#'
#' @return Vector with unique genes
#'
#' @export
getUniqueGenes <- function(x, entrez=TRUE) {

   if (!require(org.Hs.eg.db)) {
      stop("org.Hs.eg.db not installed")
   } else {
      library(org.Hs.eg.db)

      temp <-  unlist(sapply(x, function(x) strsplit(x, ";")))
      ans <- unique(temp)
      ans <- ans[!is.na(ans)]
      if (entrez)
         ans <- mapIds(org.Hs.eg.db, ans, 'ENTREZID', 'SYMBOL')
   }
   ans
}
