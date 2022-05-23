#' Exclude CpGs
#'
#' Filter 450K or EPIC methylation array with previously defined conditions (selected conditionsin config.r). In order to minimize data size we merge only the selected ethnic column
#'
#' @param cohort Dataframe with data
#' @param cpgid string with the column name that contains the CpGs identifiators.
#' @param exclude character vector with exclude parameters. Possibles values to exclude are
#' \itemize{
#'   \item{MASK_sub25_copy}{indicate whether the 25bp 3'-subsequence of the probe is non-unique}
#'   \item{MASK_sub30_copy}{indicate whether the 30bp 3'-subsequence of the probe is non-unique}
#'   \item{MASK_sub35_copy}{indicate whether the 35bp 3'-subsequence of the probe is non-unique}
#'   \item{MASK_sub40_copy}{indicate whether the 40bp 3'-subsequence of the probe is non-unique}
#'   \item{MASK_mapping}{"hether the probe is masked for mapping reason. Probes retained should have high quality (>=40 on 0-60 scale) consistent (with designed MAPINFO) mapping (for both in the case of type I) without INDELs . }
#'   \item{MASK_extBase }{Probes masked for extension base inconsistent with specified color channel (type-I) or CpG (type-II) based on mapping. }
#'   \item{MASK_typeINextBaseSwitch }{Whether the probe has a SNP in the extension base that causes a color channel switch from the official annotation (described as color-channel-switching, or CCS SNP in the reference). These probes should be processed differently than designed (by summing up both color channels instead of just the annotated color channel).}
#'   \item{MASK_snp5_common }{Whether 5bp 3'-subsequence (including extension for typeII) overlap with any of the common SNPs from dbSNP (global MAF can be under 1%). }
#'   \item{MASK_snp5_GMAF1p }{Whether 5bp 3'-subsequence (including extension for typeII) overlap with any of the SNPs with global MAF >1% . }
#'   \item{MASK_general}{ Recommended general purpose masking merged from "MASK.sub30.copy", "MASK.mapping", "MASK.extBase", "MASK.typeINextBaseSwitch" and "MASK.snp5.GMAF1p" . }
#'   \item{cpg_probes}{ cpg probes classified as "cg" in the variable named "probeType". }
#'   \item{noncpg_probes}{ non-cpg probes classified as "ch" in the variable named "probeType". }
#'   \item{control_probes}{ control probes classified as "rs" in the variable named "probeType". }
#'   \item{Unrel_450_EPIC }{Unreliable probes discordant between 450K and EPIC for blood }
#'   \item{Unrel_450_EPIC }{Unreliable probes discordant between 450K and EPIC for placenta }
#'   \item{Unrel_450_EPIC }{Unreliable probes discordant between 450K and EPIC for placenta (more restrictive) }
#'   \item{MASK_rmsk15}{ }
#'   \item{Sex}{ Keep probes targeting cpgs from sex chromosomes "chrX" and "chrY". ( CpG_chrm %in% "chrX" & CpG_chrm %in% "chrY" )}
#' }
#' @param ethnic . Ethnicity, possible values : '' (for none ethnic filter), EUR, SAS, AMR, GWD, YRI, TSI, IBS, CHS, PUR, JPT, GIH, CH_B, STU, ITU, LWK, KHV, FIN, ESN, CEU, PJL, AC_B, CLM, CDX, GBR, BE_B, PEL, MSL, MXL, ASW or GLOBAL.
#' @param artype default value '450K'. Array type Illumina 450K or Illumina EPIC arrays, possible values are '450K' or 'EPIC' respectively.
#' @param filename Optional, filename to write excluded cpgs and related information, if NULL no data are writed in a file
#' @param fileresume Optional, filename to write descriptives for removed cpgs, if NULL no data are writed in a file
#'
#' @return dataframe without the CpGs that meet the conditions
#'
#' @export
exclude_CpGs <- function(cohort, cpgid, exclude, ethnic, artype='450K', filename = NULL, fileresume = NULL)
{
   # Test if filter data exists and gets it
   if( toupper(artype) == '450K' ) {
      try(data("filter_450K") )
      filters <- filter_450K
   }else if(toupper(artype) == 'EPIC') {
      try(data("filter_EPIC"))
      filters <- filter_EPIC
   }else{
      stop( paste0( "Unknown array type ", toupper(artype) ) )
   }


   # In order to minimize data size we merge only the selected ethnic column


   # Get first ethnic column position in filters
   firstEthnicPosition <-  min(grep(paste0("MASK_snp5_.*[^GMAF1p&^common]"), colnames(filters), perl = TRUE))-1

   # Select ethnic columns not related with our data from filters (ONLY IF ethnic != '' or NA)
   if(ethnic!='' && !is.na(ethnic)) {
      # Get first ethnic column position in filters
      fieldstodelete <- grep(paste0("MASK_snp5_.*[^",ethnic,"]"),  colnames(filters)[grep(paste0("MASK_snp5_.*[^GMAF1p&^common]"), colnames(filters), perl = TRUE)], perl = TRUE) + firstEthnicPosition

   } else {
      #.Remove from code.# fieldstodelete <- grep(paste0("MASK_snp5_.*[^EUR]"),  colnames(filters)[grep(paste0("MASK_snp5_.*[^GMAF1p&^common]"), colnames(filters), perl = TRUE)], perl = TRUE) + firstEthnicPosition
      fieldstodelete <- ""
   }


   fieldstomerge <- which(!seq(1:dim(filters)[2]) %in% fieldstodelete)

   # Merge cohort with CpGs filters
   cohort <- merge(cohort, filters[,fieldstomerge], by.x= cpgid, by.y = "probeID", all.x = TRUE )


   # CpGs id to remove
   if( is.null(exclude)) {
      excludeid <- NULL
   } else {
      excludeid <- cohort[eval(parse(text=getCritera(exclude, ethnic))), cpgid]
   }


   # Report descriptive exclussions to a descriptive file
   if(!is.null(fileresume)) {
      write(sprintf('\n# %s', strrep("-",16)), file = fileresume, append = TRUE)
      write(sprintf('# Remove "problematic"  CpGs : '), file = fileresume, append = TRUE)
      write(sprintf('# %s\n', strrep("-",16)), file = fileresume, append = TRUE)
      write(sprintf('# Criteria : \n\tArray type : %s \n\tEthnia : %s \n', toupper(artype), toupper(ethnic)), file = fileresume, append = TRUE)
      write(sprintf('# Mask : %s', exclude), file = fileresume, append = TRUE)
      write(sprintf('\n'), file = fileresume, append = TRUE)
      write(sprintf('# Total CpGs in data : %d', dim(cohort)[1]), file = fileresume, append = TRUE)
      write(sprintf('# Number of excluded CpGs: %d', length(excludeid)), file = fileresume, append = TRUE)
      write(sprintf('# Total CpGs after exclusion : %d\n', (dim(cohort)[1]) - length(excludeid)), file = fileresume, append = TRUE)
      write(sprintf('# Percent excluded CpGs: %f %%\n', ((length(excludeid)/dim(cohort)[1])*100 )), file = fileresume, append = TRUE)

      suppressWarnings(
         write.table( cohort[eval(parse(text=getCritera(exclude, ethnic))),],
                      filename, col.names = TRUE, row.names = FALSE, sep = '\t', append = FALSE, dec='.'))
   }

   # # Report CpG excluded and reason to a file
   # if(!is.null(filename)) {
   #    # write(sprintf('# Criteria : %s\n', toupper(artype)), file = filename)
   #    # write(sprintf('# Total CpGs in data : %d', dim(cohort)[1]), file = filename, append = TRUE)
   #    # write(sprintf('# Number of excluded CpGs: %d', length(excludeid)), file = filename, append = TRUE)
   #    # write(sprintf('# Total CpGs after exclusion : %d\n', (dim(cohort)[1]) - length(excludeid)), file = filename, append = TRUE)
   #    # write(sprintf('# Percent excluded CpGs: %f %%\n', ((length(excludeid)/dim(cohort)[1])*100 )), file = filename, append = TRUE)
   #    # Report exclusion reason
   #    suppressWarnings(
   #       write.table( cohort[eval(parse(text=getCritera(exclude, ethnic))),],
   #                    filename, col.names = TRUE, row.names = FALSE, sep = '\t', append = FALSE, dec='.'))
   # }

   # Rmove CpGs with exclusion parameters
   cohort <- cohort[ !cohort[,cpgid] %in% excludeid,]

   warning( length(excludeid), ' CpGs have been excluded')

   return(cohort)

}




# Gets criteria dynamically from configuration variables
# function called internally by exclude_CpGs
getCritera <- function(exclude, ethnic)
{

   criter = ''
   possible_crit <- c( 'MASK_sub25_copy', 'MASK_sub30_copy', 'MASK_sub35_copy', 'MASK_sub40_copy',
                       'MASK_mapping', 'MASK_extBase', 'MASK_typeINextBaseSwitch', 'MASK_snp5_common', 'MASK_snp5_GMAF1p',
                       'MASK_general', 'cpg_probes', 'noncpg_probes', 'control_probes', 'Unrel_450_EPIC_blood', 'MASK_rmsk15',
                       'Sex', 'Unrel_450_EPIC_pla_restrict', 'Unrel_450_EPIC_pla', 'MASK_snp5_ethnic')

   if( !is.null(exclude[1]) && !is.na(exclude[1]) && exclude[1] != '')
   {
      # Test if all parameters are allowed
      if(length(which(! exclude %in% possible_crit))>=1)
         stop(paste0('Parameter(s) : ',paste(exclude[which(! exclude %in% possible_crit)], sep = ','),' not valid. Possible values are ',possible_crit ))

      # Get all exclude variables
      ##.. Works with lists -> convert list to matrix ..## exclusion.crit <- do.call(cbind, lapply( ls(patt="exclude"), get) )

      # Create formula with exclusion criteria
      ##.. Works with lists -> gets only parameters with value 'Exclude'..## criter <- paste(paste("cohort$",attributes(exclusion.crit)$dimnames[[1]][which(exclusion.crit=='Exclude')], sep = ""),"TRUE | ",sep=" == ", collapse = '')
      criter <- paste(paste("cohort$",exclude, sep = ""),"TRUE | ",sep=" == ", collapse = '')

      # If  MASK_general = 'Exclude' --> Change to :  ("MASK.sub30.copy", "MASK.mapping", "MASK.extBase", "MASK.typeINextBaseSwitch" and "MASK.snp5.GMAF1p") == TRUE
      if( length(grep("cohort$MASK_general == TRUE",criter, fixed = TRUE))>0 )
         criter <- sub("cohort$MASK_general == TRUE", " cohort$MASK_sub30_copy == TRUE | cohort$MASK_mapping == TRUE | cohort$MASK_extBase == TRUE | cohort$MASK_typeINextBaseSwitch == TRUE | cohort$MASK_snp5_GMAF1p == TRUE " , criter, fixed = TRUE)

      # If Sex = 'Exclude' --> Change to : ( CpG_chrm %in% "chrX" | CpG_chrm %in% "chrY" )
      if( length(grep("cohort$Sex == TRUE",criter, fixed = TRUE))>0 )
         criter <- sub("cohort$Sex == TRUE", " cohort$CpG_chrm %in% 'chrX' | cohort$CpG_chrm %in% 'chrY' " , criter, fixed = TRUE)

      newCpGcond <- vector()
      if( length( grep("| cohort$cpg_probes == TRUE",criter, fixed = TRUE)) >0 )
         newCpGcond <- append(newCpGcond, c('cg'))
      if( length(grep("| cohort$noncpg_probes == TRUE",criter, fixed = TRUE))>0 )
         newCpGcond <- append(newCpGcond, c('ch'))
      if( length(grep("| cohort$control_probes == TRUE",criter, fixed = TRUE))>0 )
         newCpGcond <- append(newCpGcond, c('rs'))

      if(!is.null(vector()))
      {
         # Remove invalid criteria
         criter <- sub("| cohort$cpg_probes == TRUE ", "", criter, fixed = TRUE)
         criter <- sub("| cohort$noncpg_probes == TRUE ", "", criter, fixed = TRUE)
         criter <- sub("| cohort$control_probes == TRUE ", "", criter, fixed = TRUE)
         criter <- sub("| cohort$MASK_snp5_ethnic == TRUE ", "", criter, fixed = TRUE)
         criter <- sub("cohort$cpg_probes == TRUE |", "", criter, fixed = TRUE)
         criter <- sub("cohort$noncpg_probes == TRUE |", "", criter, fixed = TRUE)
         criter <- sub("cohort$control_probes == TRUE |", "", criter, fixed = TRUE)
         criter <- sub("cohort$MASK_snp5_ethnic == TRUE |", "", criter, fixed = TRUE)


         # Add new criteria
         criter <- paste0(criter, paste0(" cohort$probeType %in% '",newCpGcond,"' | ", collapse = ''))
      }

      # # Adds ethnicity exclusion criteria
      # criter <- paste0(criter, " cohort$MASK_snp5_",ethnic," == TRUE |")

      criter <- paste0(criter, paste0(" cohort$probeType %in% '",newCpGcond,"' | ", collapse = ''))

      if( is.na(ethnic) | ethnic==''){
         # Remove extra '|' and add which function to exclusion criteria
         criter <- paste("which(",substr(criter, 0, nchar(criter)-3),")",sep="")
      }
   }

   if(!is.na(ethnic) & ethnic!=''){
      # Adds ethnicity exclusion criteria
      criter <- paste0(criter, " cohort$MASK_snp5_",ethnic," == TRUE |")

      # Remove extra '|' and add which function to exclusion criteria
      criter <- paste("which(",substr(criter, 0, nchar(criter)-2),")",sep="")
   }

   criter
}

