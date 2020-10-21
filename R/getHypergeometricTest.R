#' HyperGeometric Test for positions relative to Island
#'
#' Get HyperGeometric Test for positions relative to Island
#'
#' @param selected vector with elements that meet conditions
#' @param complete vector with all possible cases that can meet conditions
#' @param varname string. position Relative to Island name. For example : "Island", "OpenSea", "S_Shore", "N_Shore", "S_Shelf", "N_Shelf"
#'
#' @return
#'
#' @export
getHypergeometricTest <- function( significative, criteria, varname)
{
   depletion <- phyper(sum( significative == 'yes' & criteria == 'yes'), sum( criteria=='yes'),
                      sum( criteria != 'yes'), sum(significative == 'yes'), lower.tail= TRUE)

   enrichment <- phyper(sum( significative == 'yes' & criteria == 'yes') - 1, sum( criteria=='yes'),
                       sum( criteria != 'yes'), sum(significative == 'yes'), lower.tail= FALSE)

   ans <- list("depletion" = depletion,
               "enrichment" = enrichment)
   return(ans)
}
