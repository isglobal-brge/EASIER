#' Find Overlaps between two genomics ranges
#'
#' Find Overlaps between two genomics ranges and returns the coincident values. Used to search overlaps between CpGs (selected data) and Fetal Placenta (complete data)
#'
#' @param selected numerical. Variable to take in to account as significative variable, could be FDR, p-value,...
#' @param complete dataframe with all relative positions to island values to perform regression
#'
#' @return A list with the Overlap results
#' \itemize{
#'  \item{"ranges"}{A subset by overlaps}
#'  \item{"hits"}{Log-likelihood of linkage model}
#'  \item{"values"}{Coincident names}
#' }
#'
#' @export
findOverlapValues <- function(selected, complete)
{
   # Find overlaps
   ranges <- subsetByOverlaps(selected, complete)
   hits <- findOverlaps(selected, complete)

   # Find the names of states associated to each position
   idx <- subjectHits(hits)
   values <- DataFrame( States  = complete$name[idx])

   ans <- list("ranges" = ranges,
               "hits" = hits,
               "values" = values)

   return(ans)

}
