#' Restrict CpGs to array type
#'
#' Restrict CpGs to array type. In EPIC array uses only CpGs in 450K arrays
#'
#' @param CpGSites, string vector with CpG sites
#' @param artype string, Illumina array type, 450K or EPIC
#'
#' @return dataframe with data restricted to Illymina array type /450K or EPIC)
#'
#'
#' @export
restrict_CpGs_to_artype <- function(CpGSites, artype)
{

   # Get library for 450K or Epic
   if( toupper(artype) == '450K'){
      library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
      ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
   }else{
      library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
      ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
   }

   # Download data
   mergedData <-  intersect( CpGSites, rownames(ann) )
   # Merge cpgs with annotations
   # mergedData <-  merge( as.data.frame(CpGSites), ann, by.x ="CpGSites", by.y = "gene")

   return(mergedData)
}
