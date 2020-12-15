#' Get ConsensusPath Over Representation Enrichment
#'
#' Get enrichment from consensusPathDb with over representation analysis
#'
#' @param accType Dataframe with CpGs data
#' @param entityType Integer vector with column position that we are interested in
#' @param accNumbers Optional, filename to write results, if NULL no results are writed in a file
#' @param accType bolean, Optional, if true, prints text "Before QC" , else, prints "After QC"
#' @param outputfile, filename to write results, if NULL no results are writed in a file
#' @importFrom brgeEnrich cpdbOverRepresentationAnalysis
#'
#' @return list with enrichment
#'
#'
#' @export
get_consensusPdb_OverRepresentation <- function( acFSet, entityType, accNumbers, accType, outputdir = ".", outputfile = NULL)
{

   overrep <- cpdbOverRepresentationAnalysis(entityType, acFSet, accNumbers, accType)

   if(!is.null(outputfile))
   {
      if(!is.null(outputdir) & !is.na(outputdir) & outputdir!='.')
         dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
      # Output filename
      filename <- tools::file_path_sans_ext(basename(outputfile))
      # Write data to a file
      write.csv(overrep, paste0(file.path( outputdir),"/",filename,"_",acFSet,".csv"))

   }

   return( overrep)


}
