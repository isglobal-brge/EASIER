#' Significative variable OR by position Relative to Island
#'
#' Get OR by Position Relative to Island taking in to account a significative variable parameter like FDR values
#'
#' @param significative numerical. Variable to take in to account as significative variable, could be FDR, p-value,...
#' @param position vector with position Relative to Island
#' @param varname string. position Relative to Island name. For example : "Island", "OpenSea", "S_Shore", "N_Shore", "S_Shelf", "N_Shelf"
#' @filename
#'
#' @return
#'
#' @export
getRelativeIslanOR <- function(significative, position, varname)
{
   rp <- as.data.frame.matrix(table(significative, position))[,c("yes","no")]

   rp1 <- rbind(nosig = rp[1, ], sig = colSums(rp[-1, ]))[c(2,1),]
   xt <- chisq.test(rp1)
   or <- (rp1[1, 1] * rp1[2, 2]) / (rp1[1, 2] * rp1[2, 1])

   return(list("RelIsland" = varname,
               "OR" = or,
               "OR.inf" = exp(log(or)-1.96*sqrt(sum(1/rp1))),
               "OR.sup" = exp(log(or)+1.96*sqrt(sum(1/rp1))),
               "p-val" = xt$p.value))
}


