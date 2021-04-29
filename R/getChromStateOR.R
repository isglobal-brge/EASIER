#' Significative variable Regression by cromatin state
#'
#' Get regression by cromatin state taking in to account a significative variable parameter like FDR values
#'
#' @param significative numerical. Variable to take in to account as significative variable, could be FDR, p-value,...
#' @param chromstate vector. Cromatin state values
#' @param varname string. Cromatin state name. For example : "TssA","TssAFlnk","TxFlnk","TxWk","Tx","EnhG","Enh","ZNF.Rpts","Het","TssBiv","BivFlnk","EnhBiv","ReprPC","ReprPCWk" or "Quies"
#'
#' @return
#'
#' @export
getChromStateOR <- function(significative, chromstate, varname)
{

   tmp_tbl <- table(significative, chromstate)
   if(dim(tmp_tbl)[2]==1) {
      fields <- c(0,0)
      if(colnames(tmp_tbl)==TRUE){
         tmp_tbl <- cbind(tmp_tbl,fields)
      } else {
         tmp_tbl <- cbind(fields, tmp_tbl)
      }
      colnames(tmp_tbl) <- c("TRUE","FALSE")
   }

   if(dim(tmp_tbl)[1]==1) {
      fields <- c(0,0)
      if(rownames(tmp_tbl)=='no'){
         tmp_tbl <- rbind(tmp_tbl,fields)
      } else {
         tmp_tbl <- cbind(fields, tmp_tbl)
      }
      rownames(tmp_tbl) <- c("no","yes")
   }

   rp <- as.data.frame.matrix(tmp_tbl)[,c("TRUE","FALSE")]

   rp1 <- rbind(nosig = rp[1, ], sig = colSums(rp[-1, ]))[c(2,1),]
   xt <- chisq.test(rp1)
   or <- (rp1[1, 1] * rp1[2, 2]) / (rp1[1, 2] * rp1[2, 1])

   return(list("ChromStates" = varname,
               "OR" = or,
               "OR.inf" = exp(log(or)-1.96*sqrt(sum(1/rp1))),
               "OR.sup" = exp(log(or)+1.96*sqrt(sum(1/rp1))),
               "p-val" = xt$p.value))
}


