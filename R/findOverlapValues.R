#' Find Overlaps between two genomics ranges
#'
#' Find Overlaps between two genomics ranges and returns the coincident values. Used to search overlaps between CpGs (selected data) and Fetal Placenta (complete data)
#'
#' @param selected Genomic Ranges Variable that we are interested to enrich with data obtained from 'comlete' genomic ranges
#' @param complete Genomic Ranges with Fetal Placenta data 15 or 18 or other genomic ranges that overlaps with selected genomic ranges
#'
#' @return A list with the Overlap results
#' \itemize{
#'  \item{"ranges"}{A subset by overlaps}
#'  \item{"hits"}{Coincident data between selected and complete}
#'  \item{"shits"}{Subject hits}
#'  \item{"qhits"}{Query hits}
#'  \item{"values"}{Coincident names associated to each position}
#' }
#'
#' @export
findOverlapValues <- function(selected, complete,  outputdir = ".", outputfile = NULL)
{
   # Find overlaps
   ranges <- subsetByOverlaps(selected, complete)
   hits <- findOverlaps(selected, complete)

   # Find the names of states associated to each position
   idx <- subjectHits(hits)
   values <- DataFrame( States  = complete$name[idx])
   qh <- queryHits(hits)

   ans <- list("ranges" = ranges,
               "hits" = hits,
               "shits" = idx,
               "qhits" = qh,
               "values" = values)

   return(ans)

}
