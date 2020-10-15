#' Get enrichment with Molecular Signatures Database (MSigDB)
#'
#' Get enrichment with Molecular Signatures Database (MSigDB c2 - v7p1). This function gets MSigDB enrichment and save the results in  _mysmeth sufix. Data from :  Walter and Eliza Hall Institute website for local bioinformatic resources. (http://bioinf.wehi.edu.au/MSigDB/v7.1/)
#'
#' @param data array or dataframe with CpGs to perform enrichment. Dataframe option is used to filter by Bonferroni FDR or p-value, in order to filter by this fields, function needs that Bonferroni with adjustments field  to be called "Bonferroni", the field with the FDR adjustments to be called "FDR" and the field with the p-value to be called "p.value"
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
#'
#' @return A list with resulth enriched data.
#'
#'
#' @export
MSigDB_enrichment <- function( data, out, filename, artype = '450K', bn=FALSE, fdr=NA, pval=NA, all = FALSE )
{

   # Output filename
   outputfile <- tools::file_path_sans_ext(basename(filename))
   outputfile <- paste(out,outputfile,sep="/")

   if(class(data) == "character")
   {
      print("ONLY CpGs")

      sigCpGs <- as.vector(t(data))

      # Annotate CpGs
      #.FORA !!!.# annot.data <- getMappedEntrezIDs(sig.cpg=sigCpGs,  array.type = artype)
      # Hs.c2 enrichment - Only to work with a complete list of CpGs
      Hs.c2.data <- gsameth(sig.cpg=sigCpGs, collection=Hs.c2, array.type = artype)
      write.table( Hs.c2.data, paste0(outputfile,"_MSigDB.txt"), quote=F, sep="\t")

      res <- list("MSigDB" = Hs.c2.data,
                  "signif" = sigCpGs)

   } else {

      if(bn==FALSE & is.na(bn) & is.na(pval) & all==FALSE)
         stop("Nothing to enrich")

      if(! "CpGs" %in% colnames(data))
         stop("In dataframes, column with CpGs must be called CpGs")

      if(all == TRUE)
      {
         sigCpGs <- as.vector(data[, "CpGs"])

         #.FORA !!!.# annot.data <- getMappedEntrezIDs(sig.cpg=sigCpGs,  array.type = artype)
         Hs.c2.data <- gsameth(sig.cpg=sigCpGs, collection=Hs.c2, array.type = artype)
         write.table( Hs.c2.data, paste0(outputfile,"_MSigDB.txt"), quote=F, sep="\t")

         res <- list("MSigDB" = Hs.c2.data,
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
               #.FORA !!!.# annot.data <- getMappedEntrezIDs(sig.cpg=sigCpGs_bn,  all.cpg=as.vector(data[,"CpGs"]), array.type = artype)
               # MSigDB enrichment - Only to work with a complete list of CpGs
               Hs.c2.data <- gsameth(sig.cpg=sigCpGs_bn,  all.cpg=as.vector(data[,"CpGs"]), collection=Hs.c2, array.type = artype )
               write.table( Hs.c2.data, paste0(outputfile,"_MSigDB_BN.txt"), quote=F, sep="\t")

               res <- list("MSigDB.bn" = Hs.c2.data,
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
               #.FORA !!!.# annot.data <- getMappedEntrezIDs(sig.cpg=sigCpGs_fdr,  all.cpg=as.vector(data[,"CpGs"]), array.type = artype)
               # MSigDB enrichment - Only to work with a complete list of CpGs
               Hs.c2.data <- gsameth(sig.cpg=sigCpGs_fdr, all.cpg=as.vector(data[,"CpGs"]), collection=Hs.c2, array.type = artype )
               write.table( Hs.c2.data, paste0(outputfile,"_MSigDB_FDR.txt"), quote=F, sep="\t")

               if(exists("res")){
                  res[["MSigDB.fdr"]] = Hs.c2.data
                  res[["signif.fdr"]] = sigCpGs_fdr
               } else {
                  res <- list("MSigDB.fdr" = Hs.c2.data,
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
               #.FORA !!!.# annot.data <- getMappedEntrezIDs(sig.cpg=sigCpGs_pval,  all.cpg=as.vector(data[,"CpGs"]), array.type = artype)
               # MSigDB enrichment - Only to work with a complete list of CpGs
               Hs.c2.data <- gsameth(sig.cpg=sigCpGs_pval,  all.cpg=as.vector(data[,"CpGs"]), collection=Hs.c2, array.type = artype )
               write.table( Hs.c2.data, paste0(outputfile,"_MSigDB_PVAL.txt"), quote=F, sep="\t")

               if(exists("res")){
                  res[["MSigDB.pval"]] = Hs.c2.data
                  res[["signif.pval"]] = sigCpGs_pval
               } else {
                  res <- list("MSigDB.pval" = Hs.c2.data,
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
