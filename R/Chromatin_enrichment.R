#' Get enrichment with Chromatin State and classify hypo i hypermethylation states
#'
#' Get enrichment with Chromatin States Database (chrom15). This function gets Chromatin states and creates a field meth_level for hypo and Hypermethylation states if beta field is found, if not, only enrichment is performed.
#'
#' @param data array or dataframe with CpGs to perform enrichment. If not beta field found then hypo and hypermethylation field is added
#' @param out string, Path where the results should be saved
#' @param filename string, File name where the results should be saved, to this name the suffix is added depending on data to enrich
#' \itemize{
#'    \item {_chrom15}{if we are enriching a CpG vector}
#'    \item {_chrom15_BN}{if we are enriching CpGs that meet the bonferroni condition}
#'    \item {_chrom15_FDR}{if we are enriching CpGs that meet the FDR condition}
#'    \item {_chrom15_PVAL}{if we are enriching CpGs that meet the p-value condition}
#' }
#' @param cpgcol numeric or string. Column index or column name with CpGs id
#' @param bn boolean. optional. If data is a dataframe with bonferroni adjustmnet, makes enrichment with CpGs that pass Bonferroni
#' @param fdr numeric optiona. If data is a dataframe with FDR adjustmnet and fdr!=NA or NULL, makes enrichment with CpGs with fdr lower than indicated value
#' @param pval numeric optional. If data is a dataframe with p-value and pval!=NA or NULL, makes enrichment with CpGs with pval lower than indicated value
#' @param all boolean, optional. enrich all CpGs
#'
#' @return A list with resulth enriched data.
#'
#'
#' @export
Chromatin_enrichment <- function( data, out, filename, cpgcol, bn=FALSE, fdr=NA, pval=NA, all = FALSE )
{

   # Output filename
   outputfile <- tools::file_path_sans_ext(basename(filename))
   outputfile <- paste(out,outputfile,sep="/")

   if(class(data) == "character")
   {
      print("ONLY CpGs")

      sigCpGs <- as.vector(t(data))

      # Add data with chromatin15 status
      chromatin.data <- addCrom15Columns(as.data.frame(sigCpGs), 1 )
      # Recode for hiper and hypomethylation
      chromatin.data <- addCrom15Columns(as.data.frame(sigCpGs), 1 )
      write.table( chromatin.data, paste0(outputfile,"_chrom15.txt"), quote=F, sep="\t")

      res <- list("Chrom" = chromatin.data,
                  "signif" = sigCpGs)

   } else {

      if(bn==FALSE & is.na(bn) & is.na(pval) & all==FALSE)
         stop("Nothing to enrich")

      #..# if(! "CpGs" %in% colnames(data))
      #..#    stop("In dataframes, column with CpGs must be called CpGs")

      if(all == TRUE)
      {
         sigCpGs <- as.vector(data[, cpgcol])

         # Add data with chromatin15 status
         chromatin.data <- addCrom15Columns(as.data.frame(sigCpGs), "sigCpGs" )

         if( "beta" %in% colnames(data))
            chromatin.data$meth_level <- ifelse( chromatin.data$beta>=0, 'Hyper', 'Hypo')

         write.table( chromatin.data, paste0(outputfile,"_chrom15.txt"), quote=F, sep="\t")

         res <- list("Chrom" = chromatin.data,
                     "signif" = sigCpGs[,cpgcol])

      } else {

         if(bn==TRUE)
         {
            if(! "Bonferroni" %in% colnames(data))
               stop("To filter by bn the Bonferroni adjustment field must be called 'Bonferroni'")

            sigCpGs_bn <- as.vector(data[which(data$Bonferroni=='yes'), ])

            if(!is.null(sigCpGs_bn) & dim(sigCpGs_bn)[1]>0 )
            {
               # Annotate CpGs with Chromatin status
               chromatin.data <- addCrom15Columns(sigCpGs_bn, cpgcol )
               # Get methylation level
               if( "beta" %in% colnames(data))
                  chromatin.data$meth_level <- ifelse( chromatin.data$beta>=0, 'Hyper', 'Hypo')
               write.table( chromatin.data, paste0(outputfile,"_chrom15_BN.txt"), quote=F, sep="\t")

               res <- list("Chrom.bn" = chromatin.data,
                           "signif.bn" = sigCpGs_bn[,cpgcol])
            } else {
               warning("No Bonferroni significative CpGs to enrich")
            }
         }

         if(!is.na(fdr)){

            if(! "FDR" %in% colnames(data))
               stop("To filter by bn the FDR adjustment field must be called 'FDR'")

            sigCpGs_fdr <- as.vector(data[which(data$FDR < fdr), ])

            if(!is.null(sigCpGs_fdr) & dim(sigCpGs_fdr)[1]>0 )
            {
               # Chromatin Status enrichment
               chromatin.data <- addCrom15Columns(sigCpGs_fdr, cpgcol )
               # Get methylation level
               if( "beta" %in% colnames(data))
                  chromatin.data$meth_level <- ifelse( chromatin.data$beta>=0, 'Hyper', 'Hypo')
               write.table( chromatin.data, paste0(outputfile,"_chrom15_FDR.txt"), quote=F, sep="\t")

               if(exists("res")){
                  res[["Chrom.fdr"]] = chromatin.data
                  res[["signif.fdr"]] = sigCpGs_fdr[,cpgcol]
               } else {
                  res <- list("Chrom.fdr" = chromatin.data,
                              "signif.fdr" = sigCpGs_fdr[,cpgcol])
               }
            } else {
               warning("No FDR significative CpGs to enrich")
            }
         }

         if( !is.na(pval) ){

            if(! "p.value" %in% colnames(data))
               stop("To filter by bn the p-values field must be called 'p.value'")

            sigCpGs_pval <- as.vector(data[which(data$p.value < pval), ])

            if(!is.null(sigCpGs_pval) & dim(sigCpGs_pval)[1]>0 )
            {
               # Chromatin Status enrichment
               chromatin.data <- addCrom15Columns(sigCpGs_pval, cpgcol )
               # Get methylation level
               if( "beta" %in% colnames(data))
                  chromatin.data$meth_level <- ifelse( chromatin.data$beta>=0, 'Hyper', 'Hypo')
               write.table( chromatin.data, paste0(outputfile,"_chrom15_PVAL.txt"), quote=F, sep="\t")

               if(exists("res")){
                  res[["Chrom.pval"]] = chromatin.data
                  res[["signif.pval"]] = sigCpGs_pval[,cpgcol]
               } else {
                  res <- list("Chrom.pval" = chromatin.data,
                              "signif.pval" = sigCpGs_pval[,cpgcol])
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
