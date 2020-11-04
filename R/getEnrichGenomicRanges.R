#' Create GRanges from data
#'
#' Create Genomic Ranges from data
#'
#' @param seqnames string vector with seqnames
#' @param start numeric vector with start position
#' @param end numeric vector with end position
#'
#' @return Data transformed in genomic ranges with name as chr;start;end
#'
#' @export
getEnrichGenomicRanges <- function(seqnames, start, end)
{
   # Creates the name to assign to GRanges
   pmdname <- paste(seqnames, start, end, sep = ';')

   Enrich.GRange <- GRanges(
      seqnames = Rle(seqnames),
      ranges = IRanges(start, end=end),
      name = pmdname
   )
   names(Enrich.GRange) <- Enrich.GRange$name

   return(Enrich.GRange)

}
