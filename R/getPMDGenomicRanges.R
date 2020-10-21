#' Create GRanges from PMD
#'
#' Create Genomic Ranges from PMD
#'
#' @param seqnames string vector with seqnames
#' @param start numeric vector with start position
#' @param end numeric vector with end position
#'
#' @return Data transformed in genomic ranges with name as chr;start;end
#'
#' @export
getPMDGenomicRanges <- function(seqnames, start, end)
{
   # Creates the name to assign to GRanges
   pmdname <- paste(seqnames, start, end, sep = ';')

   PMD.GRange <- GRanges(
      seqnames = Rle(seqnames),
      ranges = IRanges(start, end=end),
      name = pmdname
   )
   names(PMD.GRange) <- PMD.GRange$name

   return(PMD.GRange)

}
