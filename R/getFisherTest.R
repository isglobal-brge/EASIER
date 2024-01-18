#' Get Fisher test for significative variable attending to criteria
#'
#' Get Fisher test by criteria taking in to account a significative variable parameter like FDR values
#'
#' @param significative numerical. Variable to take in to account as significative variable, could be FDR, p-value,...
#' @param criteria vector with criteria
#' @param varname string. criteria Relative to Island name. For example : "Island", "OpenSea", "S_Shore", "N_Shore", "S_Shelf", "N_Shelf" or other criteria like gene position
#'
#' @return
#'
#' @export
getFisherTest <- function(significative, criteria, varname)
{

   rp <- as.data.frame.matrix(table(significative, criteria))



   if( all(rownames(rp) %in% c('yes', 'no')) ){
      if(dim(rp)[1]<2) {
         rp <- rbind(rp,c(0,0))
         if(!'yes' %in% rownames(rp)) {
            rownames(rp) <- c('no', 'yes')
         }else {
            rownames(rp) <- c('yes', 'no')
         }
      }
   }

   if( any(rownames(rp) %in% c('Hypo-yes','Hyper-yes' )) && sum(rownames(rp) %in% c('Hypo-yes','Hyper-yes','Hypo-no','Hyper-no' ))<4 ) {
      rnames <- rownames(rp)
      rp <- rbind(rp,c(0,0))
      if(!'Hyper-yes' %in% rownames(rp)) {
         rownames(rp) <- c(rnames, 'Hyper-yes')
      }
      if(!'Hypo-yes' %in% rownames(rp)) {
         rownames(rp) <- c(rnames, 'Hypo-yes')
      }
      if(!'Hyper-no' %in% rownames(rp)) {
         rownames(rp) <- c(rnames, 'Hyper-no')
      }
      if(!'Hypo-no' %in% rownames(rp)) {
         rownames(rp) <- c(rnames, 'Hypo-no')
      }
   }


   if(dim(rp)[2]<2) {
      rp <- cbind(rp,c(0,0))
      if(!'yes' %in% colnames(rp)) {
         rp$yes <- rep(0,2)
      }else {
         rp$no <- rep(0,2)
      }
   }

   rp <- rp[,c("yes","no")]

   if(all(rownames(rp) %in% c('yes', 'no'))){
      rp1 <- rbind(nosig = rp[1, ], sig = colSums(rp[-1, ]))[c(2,1),]
      xt <- chisq.test(rp1)
      or <- (rp1[1, 1] * rp1[2, 2]) / (rp1[1, 2] * rp1[2, 1])
   } else {
      rp1 <- rp[c("Hyper-yes","Hypo-yes"), ]
      xt <- chisq.test(rp1)
      or <- (rp1[1, 1] * rp1[2, 2]) / (rp1[1, 2] * rp1[2, 1])
   }

   return(list("Data" = varname,
               "OR" = or,
               "OR.inf" = exp(log(or)-1.96*sqrt(sum(1/rp1))),
               "OR.sup" = exp(log(or)+1.96*sqrt(sum(1/rp1))),
               "p-val" = xt$p.value))
}


