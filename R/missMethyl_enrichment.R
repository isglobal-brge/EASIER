#' Get Pathways enrichment with missMethyl
#'
#' Get Pathways enrichment with missMethyl Bioconductor package. This function gets GO enrichment and KEGG enrichment from input data and save the results in a two separated files, one with GO enrichment and the other with KEGG enrichment with _GO and _KEGG sufix respectively.
#'
#' @param data array or dataframe with CpGs to perform enrichment. Dataframe option is used to filter by Bonferroni FDR or p-value, in order to filter by this fields, function needs that Bonferroni with adjustments field  to be called "Bonferroni", the field with the FDR adjustments to be called FDR and the field with the p-value to be called p.value
#' @param out string, Path where the results should be saved
#' @param filename string, File name where the results should be saved, to this name the suffix is added depending on data to enrich
#' \itemize{
#'    \item {_mysmeth}{if we are enriching a CpG vector}
#'    \item {_mysmeth_BN}{if we are enriching CpGs that meet the bonferroni condition}
#'    \item {_mysmeth_FDR}{if we are enriching CpGs that meet the FDR condition}
#'    \item {_mysmeth_PVAL}{if we are enriching CpGs that meet the p-value condition}
#' }
#' @param artype string, Illumina array type, 450K or EPIC, by default array type is 450K
#' @param bn boolean. optional. If data is a dataframe with bonferroni adjustmnet, makes enrichment with CpGs that pass Bonferroni
#' @param fdr numeric optiona. If data is a dataframe with FDR adjustmnet and fdr!=NA or NULL, makes enrichment with CpGs with fdr lower than indicated value
#' @param pval numeric optional. If data is a dataframe with p-value and pval!=NA or NULL, makes enrichment with CpGs with pval lower than indicated value
#' @param all boolean, optional. enrich all CpGs
#' @param plots boolean, optional. generate plots with GO and KEGG data
#'
#' @return File with descriptive and plots from MyssMethyl results. Results are stored in a folder MyssMethyl inside Enrichment default path
#'
#'
#' @export
missMethyl_enrichment <- function( data, out, filename, artype = '450K', bn=FALSE, fdr=NA, pval=NA, all = FALSE, outputdir = NULL, plots = FALSE )
{

   # Output filename
   outfilename <- tools::file_path_sans_ext(basename(filename))
   outdir <- paste(out,"missMethyl",sep="/")
   outfilename <- paste(outdir,outfilename,sep="/")

   print(outfilename)

   # Create dir for missMethyl files if not exists
   if(!dir.exists(outdir))
      suppressWarnings(dir.create(outdir, recursive = TRUE))

   # Data is only a CpGs vector?
   #..# if(dim(data)[1] <= 1 | dim(data)[2] <= 1)
   if(class(data) == "character")
   {
      print("ONLY CpGs")

      sigCpGs <- as.vector(t(data))

      # Annotate CpGs
      annot.data <- getMappedEntrezIDs(sig.cpg=sigCpGs,  array.type = artype)

      # GO enrichment - Only to work with a complete list of CpGs
      GO.data <- gometh(sig.cpg=sigCpGs, collection="GO", array.type = artype, plot.bias = TRUE, prior.prob = TRUE)

      # KEGG enrichment - Only to work with a complete list of CpGs
      KEGG.data <- gometh(sig.cpg=sigCpGs, collection="KEGG", array.type = artype, plot.bias = TRUE, prior.prob = TRUE)
      write_enrichment_to_file(GO.data, KEGG.data, paste0(outfilename,"_mysmeth"))

      if(plots){
         plot_missMethyl_Summary(GO.data, "ontology", "GO", "", paste0(outfilename,"_mysmeth_GO.png"));
#..#          plot_missMethyl_Summary(KEGG.data, "ontology", "KEGG", "", paste0(outfilename,"_mysmeth_KEGG.png"));
      }


      res <- list("GO" = GO.data,
                  "KEGG" = KEGG.data,
                  "signif" = sigCpGs)

   } else {

      if(bn==FALSE & is.na(bn) & is.na(pval) & all==FALSE)
         stop("Nothing to enrich")

      if(! "CpGs" %in% colnames(data))
         stop("In dataframes, column with CpGs must be called CpGs")

      if(all == TRUE)
      {
         sigCpGs <- as.vector(data[, "CpGs"])

         annot.data <- getMappedEntrezIDs(sig.cpg=sigCpGs,  array.type = artype)
         GO.data <- gometh(sig.cpg=sigCpGs, collection="GO", array.type = artype, plot.bias = TRUE, prior.prob = TRUE)
         KEGG.data <- gometh(sig.cpg=sigCpGs, collection="KEGG", array.type = artype, plot.bias = TRUE, prior.prob = TRUE)
         write_enrichment_to_file(GO.data, KEGG.data, paste0(outfilename,"_mysmeth"))
         if(plots){
            plot_missMethyl_Summary(GO.data, "ontology", "GO", "", paste0(outfilename,"_mysmeth_GO.png"));
#..#            plot_missMethyl_Summary(KEGG.data, "ontology", "KEGG", "", paste0(outfilename,"_mysmeth_KEGG.png"));
         }

         res <- list("GO" = GO.data,
                     "KEGG" = KEGG.data,
                     "signif" = sigCpGs)

      } else {

         if(bn==TRUE)
         {
            if(! "Bonferroni" %in% colnames(data))
               stop("To filter by bn the Bonferroni adjustment field must be called 'Bonferroni'")

            sigCpGs_bn <- as.vector(data[which(data$Bonferroni=='yes'), "CpGs"])

            if(!is.null(sigCpGs_bn) & length(sigCpGs_bn)>0)
            {
               # Annotate CpGs
               annot.data <- getMappedEntrezIDs(sig.cpg=sigCpGs_bn,  all.cpg=as.vector(data[,"CpGs"]), array.type = artype)

               # GO enrichment - Only to work with a complete list of CpGs
               GO.data <- gometh(sig.cpg=sigCpGs_bn,  all.cpg=as.vector(data[,"CpGs"]), collection="GO", array.type = artype,  plot.bias = TRUE, prior.prob = TRUE)
               # KEGG enrichment - Only to work with a complete list of CpGs
               KEGG.data <- gometh(sig.cpg=sigCpGs_bn, all.cpg=as.vector(data[,"CpGs"]), collection="KEGG", array.type = artype, plot.bias = TRUE, prior.prob = TRUE)

               write_enrichment_to_file(GO.data, KEGG.data, paste0(outfilename,"_mysmeth_BN"))

               if(plots){
                  plot_missMethyl_Summary(GO.data, "ontology", "GO", "", paste0(outfilename,"_mysmeth_BN_GO.png"));
                  #..#                  plot_missMethyl_Summary(KEGG.data, "ontology", "KEGG", "", paste0(outfilename,"_mysmeth_BN_KEGG.png"));
               }

               res <- list("GO.bn" = GO.data,
                           "KEGG.bn" = KEGG.data,
                           "signif.bn" = sigCpGs_bn)
            } else {
               warning("No Bonferroni significative CpGs to enrich")
            }
         }

         if(!is.na(fdr)){

            if(! "FDR" %in% colnames(data))
               stop("To filter by bn the FDR adjustment field must be called 'FDR'")

            sigCpGs_fdr <- as.vector(data[which(data$FDR < fdr), "CpGs"])

            if(!is.null(sigCpGs_fdr) & length(sigCpGs_fdr)>0) {
               # Annotate CpGs
               annot.data <- getMappedEntrezIDs(sig.cpg=sigCpGs_fdr,  all.cpg=as.vector(data[,"CpGs"]), array.type = artype)
               # GO enrichment - Only to work with a complete list of CpGs
               GO.data <- gometh(sig.cpg = sigCpGs_fdr,  all.cpg=as.vector(data[,"CpGs"]), collection="GO", array.type = artype, plot.bias = TRUE, prior.prob = TRUE)
               # KEGG enrichment - Only to work with a complete list of CpGs
               KEGG.data <- gometh(sig.cpg = sigCpGs_fdr, all.cpg=as.vector(data[,"CpGs"]), collection="KEGG", array.type = artype, plot.bias = TRUE, prior.prob = TRUE)

               write_enrichment_to_file(GO.data, KEGG.data, paste0(outfilename,"_mysmeth_FDR"))

               if(plots){
                  plot_missMethyl_Summary(GO.data, "ontology", "GO", "", paste0(outfilename,"_mysmeth_FDR_GO.png"));
                  #..#                  plot_missMethyl_Summary(KEGG.data, "ontology", "KEGG", "", paste0(outfilename,"_mysmeth_FDR_KEGG.png"));
               }

               if(exists("res")){
                  res[["GO.fdr"]] = GO.data
                  res[["KEGG.fdr"]] = KEGG.data
                  res[["signif.fdr"]] = sigCpGs_fdr
               } else {
                  res <- list("GO.fdr" = GO.data,
                              "KEGG.fdr" = KEGG.data,
                              "signif.fdr" = sigCpGs_fdr)
               }
            } else {
               warning("No FDR significative CpGs to enrich")
            }
         }

         if( !is.na(pval) ){

            if(! "p.value" %in% colnames(data))
               stop("To filter by bn the p-values field must be called 'p.value'")

            sigCpGs_pval <- as.vector(data[which(data$p.value < pval), "CpGs"])

            if( !is.null(sigCpGs_pval) & length(sigCpGs_pval)>0 ) {
               # Annotate CpGs
               annot.data <- getMappedEntrezIDs(sig.cpg=sigCpGs_pval,  all.cpg=as.vector(data[,"CpGs"]), array.type = artype)
               # GO enrichment - Only to work with a complete list of CpGs
               GO.data <- gometh(sig.cpg = sigCpGs_pval,  all.cpg=as.vector(data[,"CpGs"]), collection="GO", array.type = artype, plot.bias = TRUE, prior.prob = TRUE)
               # KEGG enrichment - Only to work with a complete list of CpGs
               KEGG.data <- gometh(sig.cpg = sigCpGs_pval, all.cpg=as.vector(data[,"CpGs"]), collection="KEGG", array.type = artype, plot.bias = TRUE, prior.prob = TRUE)

               write_enrichment_to_file(GO.data, KEGG.data, paste0(outfilename,"_mysmeth_PVAL"))

               if(plots){
                  plot_missMethyl_Summary(GO.data, "exposures", "GO", "", paste0(outfilename,"_mysmeth_pval_GO.png"));
                  #..#                  plot_missMethyl_Summary(KEGG.data, "exposures", "KEGG", "", paste0(outfilename,"_mysmeth_pval_KEGG.png"));
               }

               if(exists("res")){
                  res[["GO.pval"]] = GO.data
                  res[["KEGG.pval"]] = KEGG.data
                  res[["signif.pval"]] = sigCpGs_pval
               } else {
                  res <- list("GO.pval" = GO.data,
                              "KEGG.pval" = KEGG.data,
                              "signif.pval" = sigCpGs_pval)
               }

            }else {
               warning("No Significative CpGs by p-value to enrich")
            }
         }
      }
   }

   if(exists("res"))
      return(res)
   else
      return(NULL)
}



# Write data to file
write_enrichment_to_file <- function(GO.data, KEGG.data, outputfile)
{
   write.table( GO.data, paste0(outputfile,"_GO.txt"), quote=F, sep="\t")
   write.table( KEGG.data, paste0(outputfile,"_KEGG.txt"), quote=F, sep="\t")
}
