#' Create hapmap file for GWAMA
#'
#' Create map file for GWAMA, needed in Manhattan plots. The hapmap correspondos to Illumina data and depends on array type.
#'
#' @param artype string Illumina arrat type
#' @param outputfile string. Optional, relative or complete path with a file name for a hapmap, if not provided hapmap file is created in current directory as hapmap.map
#' @param gwama.dir string. Route to GWAMA binary
#'
#' @return Disk file with hapmap
#'
#'
#' @export
generate_hapmap_file <- function(artype, outputfile = 'hapmap.map')
{
   if( toupper(artype) == '450K' ) {
      try(data("filter_450K") )
      filters <- filter_450K
   }else {
      try(data("filter_EPIC"))
      filters <- filter_EPIC
   }

   hapmap <- as.data.frame( cbind( substr(filters$CpG_chrm,4,length(filters$CpG_chrm)),
                                   "MARKERNAME" = levels(filters$probeID),
                                   "0",
                                   filters$CpG_beg ))
   suppressWarnings( write.table( hapmap, outputfile, row.names = FALSE, col.names = FALSE , sep = '\t', append = FALSE, quote = FALSE))

}
