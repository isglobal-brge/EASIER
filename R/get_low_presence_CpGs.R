#' Get CpGs with low presence
#'
#' Returns CpGs with low presence detected after first GWAS Meta-Analysis
#'
#' @param outputfiles filename with results from GWAMA
#' @param pcentMissing per cent of missing data to be excluded
#'
#' @return array with low CpGs ids detected based on pcentMissing (variable from configuration file config_postQR)
#'
#' @export
get_low_presence_CpGs <- function(outputfiles, pcentMissing=NULL)
{
   data <- get_GWAMA_results(outputfiles[1])
   if(is.null(pcentMissing) | pcentMissing==''){
      pcentMissing = 0.8;
   }


   # Get presence by CpGs
   effectpositions <- grep("X[0-9]",names(data)) # Posicions
   dataeffects <- data[,effectpositions]
   presence <- 1-(rowSums(dataeffects == "?")/length(effectpositions))
   lowCpGs <- data[which(presence < pcentMissing),"rs_number"]

   return(lowCpGs)
}


get_GWAMA_results <- function( file )
{
   ## Readed in Remove Low Presence
   file <- paste0(file,".out")

   # Read data from gwama results
   data <- read.table(file, header = TRUE, stringsAsFactors=FALSE )

   # Split effects in columns
   lst <- strsplit(data$effects, "")
   data <- cbind(data, data.frame(matrix(unlist(lst), nrow=length(lst), byrow=T)))

   return(data)

}
