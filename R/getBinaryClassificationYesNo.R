#' Get Binary classification (yes / no)
#'
#' Get binary classificatoin from vector values taking into account a logical operator
#'
#' @param value vector with values to classify as binary "yes" / "no"
#' @param condition Condition to apply in order to classify values, "<", ">", "<=" or ">=", by default condition is "<"
#' @param limit limit value to classify data taking into account "condition" param, by default limit is 0
#'
#' @return Vector with betas length with classification in Hyper and Hypo
#'
#' @export
getBinaryClassificationYesNo <- function(value, condition = "<", limit=0)
{

   if (! condition %in% c("<", ">", "<=", ">="))
      stop("Condition must be a logical operator '<', '>', '<=' or '>=' ")

   if (!is.numeric(value)) {
      if(str_detect(value[1], ',')) {
         value <- as.numeric(sub(",", ".", value, fixed = TRUE))
      } else {
         stop("Value must be a numeric vector")
      }
   }


   # Get methylation level
   #..# binclass <- ifelse( eval(parse(text = paste(value,cond,limit))), 'yes', 'no')
   binclass <- sapply(value,function(cfdr) eval(parse(text = paste("ifelse( ",cfdr,"<", limit,", 'yes', 'no' )"))))

   return(binclass)
}
