#' Summarized States Fetal Placenta for significant, hyper and hypo in a data frame
#'
#' Summarized States Fetal Placenta for significant, hyper and hypo in a data frame
#'
#' @param significative, boolean vector with significant = 'yes' in significant values for BN, FDR, p-value ....
#' @param significativehyper boolean vector with significanthyper='yes' if CpG is significant for BN, FDR, p-value ... and CpG is hyper-methylated
#' @param significativehypo boolean vector with significanthypo='yes' if CpG is significant for BN, FDR, p-value ... and CpG is hyper-methylated
#' @param criteria string vector with values to take in to account to perform summarize, for example Chromatine States 15 Fetal Placenta
#' @param resformat string, if 'long' returns a extended table to use with plot else returns a compacted table
#' @param outputdir string. Output path to store file with results, by default results are written in current dir
#' @param outputfile string. File name to store results if no name is provided results are not written.The suffix "_States_FP" and the prefix "Porp" are added to file name provided.
#' @param plot boolean. If plot is TRUE, plot results
#'
#' @return Dataframe with summarized data (extended version)
#'
#'
#' @export
summary_HyperGeometrics_Table <- function( significative, significativehyper, significativehypo, criteria, resformat = 'long', outputdir = '.', outputfile = NULL,  plot=TRUE )
{
   crit <- unique(criteria)

   All_CpGs <-  sapply(crit, function(cr) get_proportions( ifelse(criteria == cr, "yes", "no")) )
   Sig_CpGs <-  sapply(crit, function(cr) get_conditional_proportions(significative, ifelse(criteria == cr, "yes", "no")) )
   Sig_Hyper_CpGs <-  sapply(crit, function(cr) get_conditional_proportions(significativehyper, ifelse(criteria == cr, "yes", "no")) )
   Sig_Hypo_CpGs <-  sapply(crit, function(cr) get_conditional_proportions(significativehypo, ifelse(criteria == cr, "yes", "no")) )

   tmp <- data.frame( "Criteria" = crit, "All_CpGs" = All_CpGs, "Sig_CpGs" = Sig_CpGs, "Sig_Hyper_CpGs" = Sig_Hypo_CpGs, "Sig_Hypo_CpGs" = Sig_Hyper_CpGs )

   if(resformat == 'long')
      tmp <- melt(tmp[,c("Criteria","All_CpGs","Sig_CpGs", "Sig_Hyper_CpGs", "Sig_Hypo_CpGs")],id.vars = 1)

   if(!is.null(outputfile)) {
      if(!is.null(outputdir) & !is.na(outputdir) & outputdir!='.')
         dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)

      # Output filename
      filename <- tools::file_path_sans_ext(basename(outputfile))

      write.csv(tmp, paste0(file.path( outputdir),"/Prop_",filename,"States_FP.csv"), quote=F, row.names=F, sep="\t")
   }

   if(plot){
      if(resformat != 'long'){
         dmp <- melt(tmp[,c("Criteria","All_CpGs","Sig_CpGs", "Sig_Hyper_CpGs", "Sig_Hypo_CpGs")],id.vars = 1)
         plot_ProportionHyperGeometrics( dmp, outputdir, outputfile)
      } else {
         plot_ProportionHyperGeometrics( tmp, outputdir, outputfile)
      }
   }

   return(tmp)

}


get_proportions <- function( criteria )
{
   prop <- sum( criteria == 'yes')/length(criteria)
   return(prop)
}


get_conditional_proportions <- function( significative, criteria )
{
   prop <- sum( significative == 'yes' & criteria == 'yes')/sum(significative=='yes')
   return(prop)
}
