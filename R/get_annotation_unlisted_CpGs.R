#' Get annotation from unlisted CpGs
#'
#' Get annotation (Illumina) from unlisted CpGs attending to Illumuna array type (450K or EPIC)
#'
#' @param singifCpGs vector with a list of CpGs names
#' @param artype string with artype to obtain the unsignificative CpGs
#'
#' @return
#'
#' @export
get_annotation_unlisted_CpGs <- function(singifCpGs, artype)
{
   # Get library for 450K or Epic
   if( artype == '450K'){
      library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
      ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
   }else if( artype == 'EPIC'){
      library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
      ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
   }else
      stop(paste0("Unknown ",artype))

   data_nl <-  as_tibble(ann) %>% filter( !Name %in% singifCpGs)

   return(data_nl)

}
