#' Get Entrez gene ID
#'
#' Get Entrez gene ID from gene Symbol
#'
#' @param symbols array with gene symbols
#'
#' @return Vector with Entrez gene IDs
#'
#' @export
getEntrezIdfromSymbol <- function(symbols)
{

   ### TODO : Map al symbols in eQTM database to EntrezID manually and save
   ### references in other dataset because som of then dosn't map wit org.Hs.eg.db
   ###  example : LOC148413 , CPSF3L or SLC35E2

   if (!require(org.Hs.eg.db)) {
      stop("org.Hs.eg.db not installed")
   } else {

      library(org.Hs.eg.db)

      if(is.factor(symbols)){
         warning("x must be a character, conversion applied")
         symbols <- as.character(symbols)
      }

      temp <-  unlist(sapply(symbols, function(x) strsplit(x, ";")))
      ans <- unique(temp)

      ans <- ans[!is.na(ans)]
      ans <- AnnotationDbi::select(org.Hs.eg.db,
                    keys = ans,
                    columns = c("ENTREZID", "SYMBOL"),
                    keytype = "SYMBOL")
   }

   return(ans)
}
