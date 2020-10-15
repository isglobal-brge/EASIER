#' Get Hiper and Hipo methylation data
#'
#' Get columns with Hyper and Hypo Methylation information
#'
#' @param betas vector with beta values to classify as hyper and hypo methylated
#'
#' @return Vector with betas length with classification in Hyper and Hypo
#'
#' @export
getHyperHypo <- function(betas)
{

   # Get methylation level
   meth_level <- ifelse( betas>=0, 'Hyper', 'Hypo')

   return(meth_level)
}
