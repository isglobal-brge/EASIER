#' Get eQTM enrichment
#'
#' Get Gene Ontology Enrichment analysis for the eQTM genes
#'
#' @param CpGs vector with a list of CpGs names
#' @param outputdir string. Output path to store file with results, by default results are written in current dir
#' @param outputfile string. File name to store results if no name is provided results are not written.The suffix "RegressionFDR_States" is added to file name provided.
#' @param plots boolean. If plot is TRUE, plot results
#'
#' @return gene list
#'
#' @export
geteQTMEnrichment <- function(CpGs, outputdir = ".", outputfile = NULL, plots = TRUE )
{

   # Gets all common eQTM CpGs
   eQTM_filtered <- merge(CpGs,eQTM, by.x="rs_number", by.y="CpG")
   if( nrow(eQTM_filtered) == 0 ) {
      return(list("eQTM" = NA,
                  "genes" = NA))
   }
   colnames(eQTM_filtered) <- c("rs_number", "p.value", "p.value.eQTM", "sigPair", "TC_gene")

   # get unique genes and EntrezID for this genes
   uniquevals <- getUniqueValues(eQTM_filtered$TC_gene)
   geneIds <- getEntrezIdfromSymbol(uniquevals)

   enriched <- enrichDGN(as.character(geneIds$ENTREZID))

   if(!is.null(outputfile)) {
      if(!is.null(outputdir) & !is.na(outputdir) & outputdir!='.')
         dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)

      # Output filename
      filename <- tools::file_path_sans_ext(basename(outputfile))

      write.csv(eQTM_filtered, paste0(file.path( outputdir),"/eQTM_data_",filename,".csv"))
      write.csv(uniquevals, paste0(file.path( outputdir),"/eQTM_genes_",filename,".csv"))

   }

   if(plots==TRUE)
   {

      # # Output filename
      # filename <- paste0(file.path( outputdir),"/eQTM_",tools::file_path_sans_ext(basename(outputfile)),"_barplot.pdf")
      #
      # ggsave(filename)
      #    barplot(enriched, showCategory = 30)
      # dev.off()
      #
      # # filename <- paste0(file.path( outputdir),"/eQTM_",tools::file_path_sans_ext(basename(outputfile)),"_Enrich_Map.pdf")
      # # ggsave(filename)
      # #    emapplot(enriched, pie_scale=1.5,layout="kk")
      # # dev.off()

   }

   return(list("eQTM" = eQTM_filtered,
               "genes" = geneIds))
}


# http://yulab-smu.top/clusterProfiler-book/chapter12.html#bar-plot

