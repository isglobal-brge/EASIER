#' Create hapmap file for GWAMA
#'
#' Create map file for GWAMA, needed in Manhattan plots. The hapmap correspondos to Illumina data and depends on array type.
#'
#' @param artype string Illumina arrat type
#' @param outputfile string. Optional, relative or complete path with a file name for a hapmap, if not provided hapmap file is created in current directory as hapmap.map
#' @param gwama.dir string. Route to GWAMA binary
#' @importFrom stringr str_replace
#'
#' @return Disk file with hapmap
#'
#'
#' @export
generate_hapmap_file <- function(artype, outputfile = 'hapmap.map')
{
   # browser()

   if( toupper(artype) == '450K' ) {
      library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
      map <- getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19)

   }else if( toupper(artype) == 'EPIC' ) {
      library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
      map <- getLocations(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

   }else{
      library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
      library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
      map <- unique(c(getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19),
               getLocations(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)))

   }

   hapmap <- as.data.frame( cbind( stringr::str_replace(as.vector(seqnames(map)), "chr",""),
                                   "MARKERNAME" = names(map),
                                   "0",
                                   start(map) ))

   suppressWarnings( write.table( hapmap, outputfile, row.names = FALSE, col.names = FALSE , sep = '\t', append = FALSE, quote = FALSE))

}
