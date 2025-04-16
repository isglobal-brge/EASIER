## ############ ##
##  Enrichment  ##
## ############ ##


## -------------------------------------
##  Install EASIER Package Code
## -------------------------------------
##
##  Uncomment this code to install EASIER package
#
# # Install devtools
# install.packages("devtools")
#
# # Install required packages
# devtools::source_url("https://raw.githubusercontent.com/isglobal-brge/EASIER/HEAD/installer.R")

# # Install EASIER package
# devtools::install_github("isglobal-brge/EASIER@HEAD")

##  END -  Install EASIER Package Code
## -------------------------------------


# Load package
require(EASIER)


## Develop test working directory
# setwd("~/Library/Mobile Documents/com~apple~CloudDocs/PROJECTES/Treballant/EASIER")
# setwd("/Users/mailos/tmp/proves")

########## ----------  VARIABLES DEFINED BY USER  ----------  ##########

# Set working directory to enrichment folder
setwd("<path to metaanalysis folder>/Enrichment")

# Files with CpG data to enrich may be a CpGs list or annotated GWAMA output
FilesToEnrich <- c('./GWAMA_Results/MetaA1/MetaA1_Fixed_Modif.out',
                   'GWAMA_Results/MetaA2/MetaA2_Fixed_Modif.out',
                   'toenrich/CpGstoEnrich.txt'
                   )

# Values for adjustment
BN <-  TRUE    # Use Bonferroni ?
FDR <- 0.7     # significance level for adjustment, if NA FDR is not used
pvalue <- 0.05 # significance level for p-value, if NA p-value is not used

# Array type, used : EPIC or 450K
# this data is defined for each file to analyse
artype <- c('450K', 'EPIC', 'EPIC')

# Result paths definition for QC, Meta-Analysis and Enrichment
results_folder <- 'QC_Results'
results_gwama <- '.'
results_enrich <- 'Enrichment'

# Enrichment type :  'BLOOD' or 'PLACENTA'
#     if enrichtype <- 'BLOOD' => enrichment with :
#                          Cromatine States : BLOOD (crom15)
#                          (To be implemented in future) Partially Methylated Domains (PMD) for Blood
#     if enrichtype <- 'PLACENTA' => enrichment with:
#                          Cromatine States : PLACENTA (FP_15) optionally (FP_18)
#                          Partially Methylated Domains (PMD) for Placenta
#     if enrichtype is different from 'BLOOD' and 'PLACENTA' we only get the missMethyl and MSigDB enrichment and the Unique genes list.
enrichtype <- 'PLACENTA'

IR_NCarrera <- TRUE

# Cromatine States Placenta Enrichment FP_18
# if enrichFP18 = TRUE the enrichment is performed wit FP_15 and FP_18
enrichFP18 <- FALSE

# Test to be used : 'Fisher' or 'Hypergeometric' if testdata is different no test will be performed
testdata <- 'Fisher'

# Perform eQTM enrichment
bEQTM <- TRUE

########## ----------  END VARIABLES DEFINED BY USER  ----------  ##########





## Check if we have any files to enrich and if these files exists
if (length(FilesToEnrich)>=1 & FilesToEnrich[1]!='') {
   for ( i in 1:length(FilesToEnrich))
      if (!file.exists(FilesToEnrich[i])) stop(paste0('File ',FilesToEnrich[i],' does not exsits, please check file' ))
}

## Check variables

if( ! toupper(enrichtype) %in% c('PLACENTA','BLOOD') )
   warning('Only enrichment with MyssMethyl and MSigDB will be done')

if( ! tolower(testdata) %in% c('fisher','hypergeometric') )
   warning('Wrong value for testdata variable, values must be "Fisher" or "Hypergeometric". No test will be performed ')



# Convert relative paths to absolute paths for FilesToEnrich
FilesToEnrich <- unlist(sapply(FilesToEnrich, function(file) { if(substr(file,1,1)!='.' & substr(file,1,1)!='/' & substr(file, 2, 2) != ':') file <- paste0('./',file) else file }))
FilesToEnrich <- sapply(FilesToEnrich, tools::file_path_as_absolute)

if(results_enrich!='.'){
   outputfolder <- file.path(getwd(), results_enrich )
}else{
   outputfolder <- file.path(getwd() )}


# Create dir to put results from enrichment
if(!dir.exists(outputfolder))
   suppressWarnings(dir.create(outputfolder))

setwd( outputfolder)

# Get which data we have to enrich
if (length(FilesToEnrich)>=1 & FilesToEnrich[1]!='')
{

   for (i in 1:length(FilesToEnrich)) {

      # Enrich all CpGs
      allCpGs <- FALSE

      # Get data
      data <- NULL
      data <- read.table(FilesToEnrich[i], header = TRUE, sep = "", dec = ".", stringsAsFactors = FALSE)

      # Is a CpG list only ? then read without headers and annotate data
      if(nrow(data) <= 1 || ncol(data) <= 1 || !any(c("FDR", "BN","p.value", "Bonferroni") %in% colnames(data)) ) {
         data <- read.table(FilesToEnrich[i], dec = ".") # Avoid header
         data <- as.vector(t(data))
         data <- get_annotattions(data, artype[i], FilesToEnrich[i], outputfolder )
         allCpGs <- TRUE
         data$chromosome <- substr(data$chr,4,length(data$chr))
         data$rs_number <- data$CpGs
      }else {
         if(! "rs_number" %in% colnames(data)) {
            if("CpGs" %in% colnames(data)) {
               data$rs_number = data$CpGs
            }else if("CpGId" %in% colnames(data)) {
               data$rs_number = data$CpGId
            }else {
               stop("Data must contain rs_number, CpGs or CpGId column with CpGs Ids")
            }
         }
      }

      ## -- Functional Enrichmnet
      ## ------------------------

      # Enrichment with missMethyl - GO and KEGG --> Writes results to outputfolder
      miss_enrich <- missMethyl_enrichment(data, outputfolder, FilesToEnrich[i], artype[i], BN, FDR, pvalue, allCpGs, plots = FALSE )

      # get unique genes from data
      geneUniv <- lapply( lapply(miss_enrich[grepl("signif", names(miss_enrich))], function(cpgs) { data[which(as.character(data$CpGs) %in% cpgs),]$UCSC_RefGene_Name}), getUniqueGenes)


      ## -- Online Tools

      # Enrichment with ConsensusPathDB
      #     - Consensus path http://cpdb.molgen.mpg.de/ (gene-set analysis – over-representation analysis)

      # Available FSet types :
      # 1 P     manually curated pathways from pathway databases
      # 2 N     interaction network neighborhood-based functional sets
      # 3 G2    Gene Ontology-based sets, GO level 2
      # 4 G3    Gene Ontology-based sets, GO level 3
      # 5 G4    Gene Ontology-based sets, GO level 4
      # 6 G5    Gene Ontology-based sets, GO level 5
      # 7 C     protein complex-based sets

      acFSet <- c('C', 'P', 'G2', 'G3')
      acType <- 'entrez-gene'

      # Get Enrichment
      CPDB_enrich <- lapply(names(geneUniv), function( data, accFSet, genes ) {
         print(data)
         lapply(accFSet,
                get_consensusPdb_OverRepresentation,
                entityType='genes',
                accNumbers=na.omit(as.character(eval(parse(text = paste0("genes$",data))))),
                accType=acType,
                outputdir = "ConsensusPathDB",
                outputfile = gsub(".", "_", data, fixed=TRUE) )},
         accFSet = acFSet, genes = geneUniv)

      names(CPDB_enrich) <- names(geneUniv)


      ## -- Molecular Enrichmnet
      ## -----------------------

      # Molecular Signatures Database enrichment
      msd_enrich <- MSigDB_enrichment(data, outputfolder, FilesToEnrich[i], artype[i], BN, FDR, pvalue, allCpGs)


      if("FDR" %in% colnames(data) | "Bonferroni" %in% colnames(data) |  "p.val" %in% colnames(data) | "p.value" %in% colnames(data))
      {

         ## -- Prepare data
         ## ---------------

         # Classify by Hyper and Hypo methylated
         data$meth_state <- getHyperHypo(data$beta) # Classify methylation into Hyper and Hypo


         if("FDR" %in% colnames(data) & !is.na(FDR) )
         {
            # Add column bFDR to data for that CpGs that accomplish with FDR
            data$bFDR <- getBinaryClassificationYesNo(data$FDR, "<", FDR) # Classify fdr into "yes" and no taking into account FDR significance level

            # CpGs FDR and Hyper and Hypo respectively
            FDR_Hyper <- ifelse(data$bFDR == 'yes' & data$meth_state=='Hyper', "yes", "no")
            FDR_Hypo <- ifelse(data$bFDR == 'yes' & data$meth_state=='Hypo', "yes", "no")
            FDR_Hyper_Hypo <- ifelse(data$bFDR == 'yes' & data$meth_state=='Hypo', "Hypo-yes",
                                      ifelse(data$bFDR == 'no' & data$meth_state=='Hypo', "Hypo-no",
                                             ifelse(data$bFDR == 'yes' & data$meth_state=='Hyper', "Hyper-yes", "Hyper-no" ) ) )
         }

         if("Bonferroni" %in% colnames(data) & BN==TRUE)
         {
            # CpGs Bonferroni and Hyper and Hypo respectively
            BN_Hyper <- ifelse(data$Bonferroni == 'yes' & data$meth_state=='Hyper', "yes", "no")
            BN_Hypo <- ifelse(data$Bonferroni == 'yes' & data$meth_state=='Hypo', "yes", "no")
            BN_Hyper_Hypo <- ifelse(data$Bonferroni == 'yes' & data$meth_state=='Hypo', "Hypo-yes",
                                    ifelse(data$Bonferroni == 'no' & data$meth_state=='Hypo', "Hypo-no",
                                           ifelse(data$Bonferroni == 'yes' & data$meth_state=='Hyper', "Hyper-yes", "Hyper-no" ) ) )
         }


         if(("p.val" %in% colnames(data) | "p.value" %in% colnames(data)) & !is.na(pvalue) )
         {

            # Add column bpval to data for that CpGs that accomplish with FDR
            data$bpval <- getBinaryClassificationYesNo(data$p.value, "<", pvalue) # Classify fdr into "yes" and no taking into account FDR significance level

            # CpGs FDR and Hyper and Hypo respectively
            pval_Hyper <- ifelse(data$bpval == 'yes' & data$meth_state=='Hyper', "yes", "no")
            pval_Hypo <- ifelse(data$bpval == 'yes' & data$meth_state=='Hypo', "yes", "no")
            pval_Hyper_Hypo <- ifelse(data$bpval == 'yes' & data$meth_state=='Hypo', "Hypo-yes",
                                      ifelse(data$bpval == 'no' & data$meth_state=='Hypo', "Hypo-no",
                                             ifelse(data$bpval == 'yes' & data$meth_state=='Hyper', "Hyper-yes", "Hyper-no" ) ) )
         }


         ## TODO: Simplify this code with only one function x option (fisher - Geometric // BN - FDR )

         # For FDR
         if("FDR" %in% colnames(data) & !is.na(FDR) )
         {

            ## --  CpG Gene position
            ## ---------------------

            # Get descriptives
            get_descriptives_GenePosition(data$UCSC_RefGene_Group, data$bFDR , "FDR", outputdir = "GenePosition/Fisher_FDR_Desc", outputfile = FilesToEnrich[i])

            if( tolower(testdata) =='fisher') {
               ## --  Fisher Test - Gene position - FDR, FDR_hyper and FDR_hypo
               GenePosition <- getAllFisherTest(data$bFDR, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_FDR", outputfile = FilesToEnrich[i], plots = TRUE )
               GenePosition_hyper <- getAllFisherTest(FDR_Hyper, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_FDRHyper", outputfile = FilesToEnrich[i], plots = TRUE )
               GenePosition_hypo <- getAllFisherTest(FDR_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_FDRHypo", outputfile = FilesToEnrich[i], plots = TRUE )
               GenePosition_HyperHypo <- getAllFisherTest(FDR_Hyper_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_FDRHyperHypo", outputfile = FilesToEnrich[i], plots = TRUE )
            }
            else if ( tolower(testdata) =='hypergeometric') {
               ## --  HyperGeometric Test - Island relative position - FDR, FDR_hyper and FDR_hypo (for Depletion and Enrichment)
               GenePosition <- getAllHypergeometricTest(data$bFDR, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_FDR", outputfile = FilesToEnrich[i])
               GenePosition_hyper <- getAllHypergeometricTest(FDR_Hyper, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_FDRHyper", outputfile = FilesToEnrich[i])
               GenePosition_hypo <- getAllHypergeometricTest(FDR_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_FDRHypo", outputfile = FilesToEnrich[i])
            }

            plot_TestResults_Collapsed(list(fdr = GenePosition, fdr_hypo = GenePosition_hypo, fdr_hyper = GenePosition_hyper),
                                       outputdir = "GenePosition", outputfile = FilesToEnrich[i], main = )

            ## --  CpG Island relative position
            ## --------------------------------

            # Get descriptives
            get_descriptives_RelativetoIsland(data$Relation_to_Island, data$bFDR , "FDR", outputdir = "RelativeToIsland/Fisher_FDR_RelativeToIsland", outputfile = FilesToEnrich[i])

            if( tolower(testdata) =='fisher') {
               ## --  Fisher Test - Position Relative to Island - FDR, FDR_hyper and FDR_hypo
               relative_island <- getAllFisherTest(data$bFDR, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_FDR", outputfile = FilesToEnrich[i], plots = TRUE )
               relative_island_hyper <- getAllFisherTest(FDR_Hyper, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_FDRHyper", outputfile = FilesToEnrich[i], plots = TRUE )
               relative_island_hypo <- getAllFisherTest(FDR_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_FDRHypo", outputfile = FilesToEnrich[i], plots = TRUE )
               relative_island_hyperhypo <- getAllFisherTest(FDR_Hyper_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_FDRHyperHypo", outputfile = FilesToEnrich[i], plots = TRUE )
            }
            else if ( tolower(testdata) =='hypergeometric') {
               ## --  HyperGeometric Test - Gene position - FDR, FDR_hyper and FDR_hypo (for Depletion and Enrichment)
               relative_island <- getAllHypergeometricTest(data$bFDR, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_FDR", outputfile = FilesToEnrich[i])
               relative_island_hyper <- getAllHypergeometricTest(FDR_Hyper, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_FDRHyper", outputfile = FilesToEnrich[i])
               relative_island_hypo <- getAllHypergeometricTest(FDR_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_FDRHypo", outputfile = FilesToEnrich[i])
            }
            plot_TestResults_Collapsed(list(bn = relative_island, bn_hypo = relative_island_hypo, bn_hyper = relative_island_hyper),
                                       outputdir = "RelativeToIsland", outputfile = FilesToEnrich[i], main = )


         }

         # For Bonferroni
         if("Bonferroni" %in% colnames(data) & BN==TRUE)
         {

            ## --  CpG Gene position
            ## ---------------------

            get_descriptives_GenePosition(data$UCSC_RefGene_Group, data$Bonferroni, "Bonferroni", outputdir = "GenePosition/Fisher_BN_Desc", outputfile = FilesToEnrich[i])

            # For BN
            if( tolower(testdata) =='fisher') {
               ## --  Fisher Test - Gene position - FDR, FDR_hyper and FDR_hypo
               GenePosition <- getAllFisherTest(data$Bonferroni, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_BN", outputfile = FilesToEnrich[i], plots = TRUE )
               GenePosition_hyper <- getAllFisherTest(BN_Hyper, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_BNHyper", outputfile = FilesToEnrich[i], plots = TRUE )
               GenePosition_hypo <- getAllFisherTest(BN_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_BNHypo", outputfile = FilesToEnrich[i], plots = TRUE )
               GenePosition_hyperhypo <- getAllFisherTest(BN_Hyper_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_BNHyperHypo", outputfile = FilesToEnrich[i], plots = TRUE )
            }
            else if ( tolower(testdata) =='hypergeometric') {
               ## --  HyperGeometric Test - Island relative position - FDR, FDR_hyper and FDR_hypo (for Depletion and Enrichment)
               GenePosition <- getAllHypergeometricTest(data$Bonferroni, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_BN", outputfile = FilesToEnrich[i])
               GenePosition_hyper <- getAllHypergeometricTest(BN_Hyper, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_BNHyper", outputfile = FilesToEnrich[i])
               GenePosition_hypo <- getAllHypergeometricTest(BN_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_BNHypo", outputfile = FilesToEnrich[i])
            }

            plot_TestResults_Collapsed(list(fdr = GenePosition, fdr_hypo = GenePosition_hypo, fdr_hyper = GenePosition_hyper),
                                       outputdir = "GenePosition", outputfile = FilesToEnrich[i], main = )

            ## --  CpG Island relative position
            ## --------------------------------

            # Get descriptives
            get_descriptives_RelativetoIsland(data$Relation_to_Island, data$Bonferroni, "Bonferroni", outputdir = "RelativeToIsland/Fisher_BN_RelativeToIsland", outputfile = FilesToEnrich[i])

            if( tolower(testdata) =='fisher') {
               ## --  Fisher Test - Position Relative to Island - FDR, FDR_hyper and FDR_hypo
               relative_island <- getAllFisherTest(data$Bonferroni, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_BN", outputfile = FilesToEnrich[i], plots = TRUE )
               relative_island_hyper <- getAllFisherTest(BN_Hyper, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_BNHyper", outputfile = FilesToEnrich[i], plots = TRUE )
               relative_island_hypo <- getAllFisherTest(BN_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_BNHypo", outputfile = FilesToEnrich[i], plots = TRUE )
               relative_island_hyperhypo <- getAllFisherTest(BN_Hyper_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_BNHyperHypo", outputfile = FilesToEnrich[i], plots = TRUE )
            }
            else if ( tolower(testdata) =='hypergeometric') {
               ## --  HyperGeometric Test - Gene position - FDR, FDR_hyper and FDR_hypo (for Depletion and Enrichment)
               relative_island <- getAllHypergeometricTest(data$Bonferroni, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_BN", outputfile = FilesToEnrich[i])
               relative_island_hyper <- getAllHypergeometricTest(BN_Hyper, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_BNHyper", outputfile = FilesToEnrich[i])
               relative_island_hypo <- getAllHypergeometricTest(BN_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_BNHypo", outputfile = FilesToEnrich[i])
            }
            plot_TestResults_Collapsed(list(bn = relative_island, bn_hypo = relative_island_hypo, bn_hyper = relative_island_hyper),
                                       outputdir = "RelativeToIsland", outputfile = FilesToEnrich[i], main = )

         }

         # For pvalue
         if( ("p.val" %in% colnames(data) | "p.value" %in% colnames(data)) & !is.na(pvalue) )
         {

            ## --  CpG Gene position
            ## ---------------------

            # Get descriptives
            get_descriptives_GenePosition(data$UCSC_RefGene_Group, data$bpval , "p.value", outputdir = "GenePosition/Fisher_pval_Desc", outputfile = FilesToEnrich[i])

            if( tolower(testdata) =='fisher') {
               ## --  Fisher Test - Gene position - pval, pval_hyper and pval_hypo
               GenePosition <- getAllFisherTest(data$bpval, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_pval", outputfile = FilesToEnrich[i], plots = TRUE )
               GenePosition_hyper <- getAllFisherTest(pval_Hyper, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_pvalHyper", outputfile = FilesToEnrich[i], plots = TRUE )
               GenePosition_hypo <- getAllFisherTest(pval_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_pvalHypo", outputfile = FilesToEnrich[i], plots = TRUE )
               GenePosition_hyperhypo <- getAllFisherTest(pval_Hyper_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_pvalHyperHypo", outputfile = FilesToEnrich[i], plots = TRUE )
            }
            else if ( tolower(testdata) =='hypergeometric') {
               ## --  HyperGeometric Test - Island relative position - pval, pval_hyper and pval_hypo (for Depletion and Enrichment)
               GenePosition <- getAllHypergeometricTest(data$bpval, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_pval", outputfile = FilesToEnrich[i])
               GenePosition_hyper <- getAllHypergeometricTest(pval_Hyper, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_pvalHyper", outputfile = FilesToEnrich[i])
               GenePosition_hypo <- getAllHypergeometricTest(pval_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_pvalHypo", outputfile = FilesToEnrich[i])
            }

            plot_TestResults_Collapsed(list(pval = GenePosition, pval_hypo = GenePosition_hypo, pval_hyper = GenePosition_hyper),
                                       outputdir = "GenePosition", outputfile = FilesToEnrich[i], main = )

            ## --  CpG Island relative position
            ## --------------------------------

            # Get descriptives
            get_descriptives_RelativetoIsland(data$Relation_to_Island, data$bpval , "p.value", outputdir = "RelativeToIsland/Fisher_pval_RelativeToIsland", outputfile = FilesToEnrich[i])

            if( tolower(testdata) =='fisher') {
               ## --  Fisher Test - Position Relative to Island - pval, pval_hyper and pval_hypo
               relative_island <- getAllFisherTest(data$bpval, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_pval", outputfile = FilesToEnrich[i], plots = TRUE )
               relative_island_hyper <- getAllFisherTest(pval_Hyper, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_pvalHyper", outputfile = FilesToEnrich[i], plots = TRUE )
               relative_island_hypo <- getAllFisherTest(pval_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_pvalHypo", outputfile = FilesToEnrich[i], plots = TRUE )
               relative_island_hyperhypo <- getAllFisherTest(pval_Hyper_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_pvalHyperHypo", outputfile = FilesToEnrich[i], plots = TRUE )
            }
            else if ( tolower(testdata) =='hypergeometric') {
               ## --  HyperGeometric Test - Gene position - pval, pval_hyper and pval_hypo (for Depletion and Enrichment)
               relative_island <- getAllHypergeometricTest(data$bpval, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_pval", outputfile = FilesToEnrich[i])
               relative_island_hyper <- getAllHypergeometricTest(pval_Hyper, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_pvalHyper", outputfile = FilesToEnrich[i])
               relative_island_hypo <- getAllHypergeometricTest(pval_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_pvalHypo", outputfile = FilesToEnrich[i])
            }
            plot_TestResults_Collapsed(list(bn = relative_island, bn_hypo = relative_island_hypo, bn_hyper = relative_island_hyper),
                                       outputdir = "RelativeToIsland", outputfile = FilesToEnrich[i], main = )


         }


      } else {

         ## -- Prepare data
         ## ---------------

         # * Create dataframe with non significative CpGs attending to artype ('EPIC' or '450K')
         #     - Get gene position tests
         #     - Get CpG Island relative position ests

         # Get
         unsignif_df <- get_annotation_unlisted_CpGs(data$rs_number, artype[i])

         data$signif <- 'yes'
         unsignif_df$signif <- 'no'

         data <- rbind(data[,which(colnames(data) %in% colnames(as.data.frame(unsignif_df)))], as.data.frame(unsignif_df) )
         data$rs_number <- data$Name


         ## TODO: Simplify this code with only one function x option (fisher - Geometric // BN - FDR )

         ## --  CpG Gene position
         ## ---------------------

         # Get descriptives
         get_descriptives_GenePosition(data$UCSC_RefGene_Group, data$signif , "CpGlist", outputdir = "GenePosition/Fisher_CpGlist_Desc", outputfile = FilesToEnrich[i])

         if( tolower(testdata) =='fisher') {
            GenePosition <- getAllFisherTest(data$signif, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_CpGlist", outputfile = FilesToEnrich[i], plots = TRUE )
         }else if ( tolower(testdata) =='hypergeometric') {
            GenePosition <- getAllHypergeometricTest(data$signif, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_CpGlist", outputfile = FilesToEnrich[i])
         }

         #..# plot_RelativetoIsland(GenePosition, outputdir = "GenePosition", outputfile = paste0("CpGlist_",FilesToEnrich[i]), main = )
         plot_GenePosition(GenePosition, outputdir = "GenePosition", outputfile = paste0("CpGlist_",FilesToEnrich[i]), main = )

         ## --  CpG Island relative position
         ## --------------------------------

         # Get descriptives
         get_descriptives_RelativetoIsland(data$Relation_to_Island, data$signif , "CpGlist", outputdir = "RelativeToIsland/Fisher_CpGlist_RelativeToIsland", outputfile = FilesToEnrich[i])

         if( tolower(testdata) =='fisher') {
            relative_island <- getAllFisherTest(data$signif, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_CpGlist", outputfile = FilesToEnrich[i], plots = TRUE )
         } else {
            relative_island <- getAllHypergeometricTest(data$signif, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_CpGlist", outputfile = FilesToEnrich[i])
         }

         plot_OR(relative_island, outputdir = "RelativeToIsland", outputfile = paste0("CpGlist_",FilesToEnrich[i]), main = )

         #..# plot_TestResults_Collapsed(list(relat = relative_island),
         #..#                            outputdir = "RelativeToIsland", outputfile = paste0("CpGlist_",FilesToEnrich[i]), main = )

      }





      ## ----------------------
      ## -- BLOOD ENRICHMENT --
      ## ----------------------

      if ( toupper(enrichtype) == 'BLOOD' )
      {
         ## --  ROADMAP  -  Metilation in Cromatine States - BLOOD
         ## -------------------------------------------------------
         ##       Analysis of methylation changes in the different chromatin states (CpGs are diff meth in some states and others don't)

         # Prepare data
         crom_data <- addCrom15Columns(data, "rs_number") # Adds chromatine state columns

         # Columns with chromatin status information :
         ChrStatCols <- c("TssA","TssAFlnk","TxFlnk","TxWk","Tx","EnhG","Enh","ZNF.Rpts","Het","TssBiv","BivFlnk","EnhBiv","ReprPC","ReprPCWk","Quies")

         if("FDR" %in% colnames(data) | "Bonferroni" %in% colnames(data) | "p.val" %in% colnames(data) & (BN==TRUE | !is.na(FDR) | !is.na(pvalue)))
         {

            if( !is.na(FDR) ) {
               chrom_states_fdr <- getAllChromStateOR( crom_data$bFDR, crom_data[,ChrStatCols], outputdir = "ChrSates_15_Blood/OR_FDR", outputfile = FilesToEnrich[i], plots = TRUE )
               chrom_states_fdr_hyper <- getAllChromStateOR( FDR_Hyper[which(data$rs_number %in% crom_data$rs_number)], crom_data[,ChrStatCols], outputdir = "ChrSates_15_Blood/OR_FDRHyper", outputfile = FilesToEnrich[i], plots = TRUE )
               chrom_states_fdr_hypo <- getAllChromStateOR( FDR_Hypo[which(data$rs_number %in% crom_data$rs_number)], crom_data[,ChrStatCols], outputdir = "ChrSates_15_Blood/OR_FDRHypo", outputfile = FilesToEnrich[i], plots = TRUE )

               plot_TestResults_Collapsed(list(fdr = chrom_states_fdr, fdr_hypo = chrom_states_fdr_hypo, fdr_hyper = chrom_states_fdr_hyper),
                                          outputdir = "ChrSates_15_Blood", outputfile = FilesToEnrich[i], main = )
            }
            if ( BN == TRUE) {
               chrom_states_bn <- getAllChromStateOR( crom_data$Bonferroni, crom_data[,ChrStatCols], outputdir = "ChrSates_15_Blood/OR_BN", outputfile = FilesToEnrich[i], plots = TRUE )
               chrom_states_bn_hyper <- getAllChromStateOR( BN_Hyper[which(data$rs_number %in% crom_data$rs_number)], crom_data[,ChrStatCols], outputdir = "ChrSates_15_Blood/OR_BNHyper", outputfile = FilesToEnrich[i], plots = TRUE )
               chrom_states_bn_hypo <- getAllChromStateOR( BN_Hypo[which(data$rs_number %in% crom_data$rs_number)], crom_data[,ChrStatCols], outputdir = "ChrSates_15_Blood/OR_BNHypo", outputfile = FilesToEnrich[i], plots = TRUE )

               plot_TestResults_Collapsed(list(bn = chrom_states_bn, bn_hypo = chrom_states_bn_hypo, bn_hyper = chrom_states_bn_hyper),
                                          outputdir = "ChrSates_15_Blood", outputfile = FilesToEnrich[i], main = )
            }
            if( !is.na(pvalue) ) {
               chrom_states_pval <- getAllChromStateOR( crom_data$bpval, crom_data[,ChrStatCols], outputdir = "ChrSates_15_Blood/OR_pval", outputfile = FilesToEnrich[i], plots = TRUE )
               chrom_states_pval_hyper <- getAllChromStateOR( pval_Hyper[which(data$rs_number %in% crom_data$rs_number)], crom_data[,ChrStatCols], outputdir = "ChrSates_15_Blood/OR_pvalHyper", outputfile = FilesToEnrich[i], plots = TRUE )
               chrom_states_pval_hypo <- getAllChromStateOR( pval_Hypo[which(data$rs_number %in% crom_data$rs_number)], crom_data[,ChrStatCols], outputdir = "ChrSates_15_Blood/OR_pvalHypo", outputfile = FilesToEnrich[i], plots = TRUE )

               plot_TestResults_Collapsed(list(pval = chrom_states_pval, pval_hypo = chrom_states_pval_hypo, pval_hyper = chrom_states_pval_hyper),
                                          outputdir = "ChrSates_15_Blood", outputfile = FilesToEnrich[i], main = )
            }
         } else {

            chrom_states <- getAllChromStateOR( crom_data$signif, crom_data[,ChrStatCols], outputdir = "ChrSates_15_Blood/OR_CpGlist", outputfile = FilesToEnrich[i], plots = TRUE )
            #..# plot_chromosomestate(chrom_states, outputdir = "ChrSates_15_Blood", outputfile = FilesToEnrich[i], main = )
         }


         # Add bEQTM enrichment
         if(bEQTM == TRUE)
         {
            if( !is.na(FDR) ) {
               geteQTMEnrichment(crom_data[which(crom_data$bFDR=='yes'),c("rs_number", "p.value")], outputdir = "eQTM_Blood/FDR", outputfile = FilesToEnrich[i], plots = TRUE )
            }
            if ( BN == TRUE) {
               geteQTMEnrichment(crom_data[which(crom_data$Bonferroni =='yes'),c("rs_number", "p.value")], outputdir = "eQTM_Blood/BN", outputfile = FilesToEnrich[i], plots = TRUE )
            }
            if( !is.na(pvalue) ) {
               geteQTMEnrichment(crom_data[which(crom_data$bpval =='yes'),c("rs_number", "p.value")], outputdir = "eQTM_Blood/pval", outputfile = FilesToEnrich[i], plots = TRUE )
            }
         }



      }


      ## ------------------------------------------------
      ## -- IMPRINTING REGIONS - NATALIA CARRERA PAPER --
      ## ------------------------------------------------

      if ( IR_NCarrera == TRUE )
      {

         ## -- Imprinting Regions CARRERA .
         ## ------------------------------------------------

         # Adds rs_number column if not in dataframe
         if(! "rs_number" %in% colnames(data)){
            data$rs_number <- data$CpGs}

         # Convert to Genomic Ranges
         data.GRange <- GRanges(
            seqnames = Rle(data$chr),
            ranges=IRanges(data$pos, end=data$pos),
            name=data$rs_number,
            chr=data$chr,
            pos=data$pos
         )
         names(data.GRange) <- data.GRange$name

         # Create genomic ranges from IRC data
         IRC.GRange <- getEnrichGenomicRanges(IRCarrera_Ev1$ICR_chr, IRCarrera_Ev1$ICR_start, IRCarrera_Ev1$ICR_stop)

         # Find overlaps between CpGs and IRC (find subject hits, query hits )
         overIRC <- findOverlapValues(data.GRange, IRC.GRange )

         #Create a data.frame with CpGs and IRCs information
         mdata <- cbind.data.frame(DataFrame(CpG = data.GRange$name[overIRC$qhits]), DataFrame(IRC = IRC.GRange$name[overIRC$shits]))

         # Merge with results from meta-analysis (A2)
         crom_data <- merge(crom_data, mdata, by.x="rs_number", by.y="CpG",all=T)

         # CpGs with IRC as NA
         IRC_NaN <- ifelse(is.na(crom_data$IRC.y),'IsNA','NotNA' )

         if("FDR" %in% colnames(data) | "Bonferroni" %in% colnames(data) & (BN==TRUE | !is.na(FDR) ))
         {

            if( tolower(testdata) =='fisher') {

               if( !is.na(FDR) ) {
                  ## --  HyperGeometric Test - IRC - FDR,  FDR_hyper and FDR_hypo  (Full data ) (Depletion and Enrichment)
                  IRC_fdr <- getAllFisherTest(crom_data$bFDR, IRC_NaN, outputdir = "IR_Carrera/Fisher_FDR", outputfile = FilesToEnrich[i])
                  IRC_fdrhyper <- getAllFisherTest(FDR_Hyper, IRC_NaN, outputdir = "IR_Carrera/Fisher_FDRHyper", outputfile = FilesToEnrich[i])
                  IRC_fdrhypo <- getAllFisherTest(FDR_Hypo, IRC_NaN, outputdir = "IR_Carrera/Fisher_FDRHypo", outputfile = FilesToEnrich[i])
               }
               if ( BN == TRUE) {
                  ## --  HyperGeometric Test - IRC - BN,  BN_hyper and BN_hypo  (Full data ) (Depletion and Enrichment)
                  IRC_bn <- getAllFisherTest(crom_data$Bonferroni, IRC_NaN, outputdir = "IR_Carrera/Fisher_BN", outputfile = FilesToEnrich[i])
                  IRC_bnhyper <- getAllFisherTest(BN_Hyper, IRC_NaN, outputdir = "IR_Carrera/Fisher_BNHyper", outputfile = FilesToEnrich[i])
                  IRC_bnhypo <- getAllFisherTest(BN_Hypo, IRC_NaN, outputdir = "IR_Carrera/Fisher_BNHypo", outputfile = FilesToEnrich[i])
               }

            } else if ( tolower(testdata) =='hypergeometric') {

               if( !is.na(FDR) ) {
                  ## --  HyperGeometric Test - IRC - FDR,  FDR_hyper and FDR_hypo  (Full data ) (Depletion and Enrichment)
                  IRC_fdr <- getAllHypergeometricTest(crom_data$bFDR, IRC_NaN, outputdir = "IR_Carrera/HyperG_FDR", outputfile = FilesToEnrich[i])
                  IRC_fdrhyper <- getAllHypergeometricTest(FDR_Hyper, IRC_NaN, outputdir = "IR_Carrera/HyperG_FDRHyper", outputfile = FilesToEnrich[i])
                  IRC_fdrhypo <- getAllHypergeometricTest(FDR_Hypo, IRC_NaN, outputdir = "IR_Carrera/HyperG_FDRHypo", outputfile = FilesToEnrich[i])
                  # Summary
                  resdata <- summary_HyperGeometrics_Table( crom_data$bFDR, FDR_Hyper, FDR_Hypo, IRC_NaN, outputdir = "IR_Carrera/Summary_HyperG_FDR", outputfile = FilesToEnrich[i], plot = TRUE )
               }

               if ( BN == TRUE) {
                  ## --  HyperGeometric Test - IRC - BN,  BN_hyper and BN_hypo  (Full data ) (Depletion and Enrichment)
                  IRC_bn <- getAllHypergeometricTest(crom_data$Bonferroni, IRC_NaN, outputdir = "IR_Carrera/HyperG_BN", outputfile = FilesToEnrich[i])
                  IRC_bnhyper <- getAllHypergeometricTest(BN_Hyper, IRC_NaN, outputdir = "IR_Carrera/HyperG_BNHyper", outputfile = FilesToEnrich[i])
                  IRC_bnhypo <- getAllHypergeometricTest(BN_Hypo, IRC_NaN, outputdir = "IR_Carrera/HyperG_BNHypo", outputfile = FilesToEnrich[i])
                  # Summary
                  resdata <- summary_HyperGeometrics_Table( crom_data$Bonferroni, BN_Hyper, BN_Hypo, IRC_NaN, outputdir = "IR_Carrera/Summary_HyperG_BN", outputfile = FilesToEnrich[i], plot = TRUE )
               }
            }

            if ( !is.na(FDR) )  {
               ## --  Plot collapsed data HyperGeometric Test - States15_FP - FDR
               plot_TestResults_Collapsed(list(fdr = IRC_fdr, fdr_hypo = IRC_fdrhypo, fdr_hyper = IRC_fdrhyper),
                                          outputdir = "IR_Carrera", outputfile = FilesToEnrich[i])
            }

            if ( BN == TRUE) {
               ## --  Plot collapsed data HyperGeometric Test - States15_FP - BN
               plot_TestResults_Collapsed(list(bn = IRC_bn, bn_hypo = IRC_bnhypo, bn_hyper = IRC_bnhyper),
                                          outputdir = "IR_Carrera", outputfile = FilesToEnrich[i])
            }

         } else {

            if( tolower(testdata) =='fisher') {
               IRC <- getAllFisherTest(crom_data$signif, IRC_NaN, outputdir = "IR_Carrera/Fisher_CpGlist", outputfile = FilesToEnrich[i])
            } else if ( tolower(testdata) =='hypergeometric') {
               IRC <- getAllHypergeometricTest(crom_data$signif, IRC_NaN, outputdir = "IR_Carrera/HyperG_CpGlist", outputfile = FilesToEnrich[i])
            }
         }

      }



      ## -------------------------
      ## -- PLACENTA ENRICHMENT --
      ## -------------------------

      if ( toupper(enrichtype) == 'PLACENTA' )
      {

         ## -- ROADMAP  -  Regulatory feature enrichment analysis - PLACENTA
         ## -----------------------------------------------------------------

         # Adds rs_number column if not in dataframe
         if(! "rs_number" %in% colnames(data)){
            data$rs_number <- data$CpGs}

         # Convert to Genomic Ranges
         data.GRange <- GRanges(
            seqnames = Rle(data$chr),
            ranges=IRanges(data$pos, end=data$pos),
            name=data$rs_number,
            chr=data$chr,
            pos=data$pos
         )
         names(data.GRange) <- data.GRange$name

         # Find overlaps between CpGs and Fetal Placenta (States 15 and 18)
         over15 <- findOverlapValues(data.GRange, FP_15_E091 )

         if (enrichFP18 == TRUE){
            over18 <- findOverlapValues(data.GRange, FP_18_E091 )
            # Add states 15 and 18 to data.GRange file and write to a file : CpGs, state15 and state18
            data.chrstates <- c(mcols(over15$ranges), over15$values, over18$values)
            colnames(data.chrstates)[grep("States",colnames(data.chrstates))] <-  c("States15_FP", "States18_FP")
         } else {
            # Add states 15 to data.GRange file and write to a file : CpGs, state15
            data.chrstates <- c(mcols(over15$ranges), over15$values)
            colnames(data.chrstates)[grep("States",colnames(data.chrstates))] <-  c("States15_FP")
         }

         # Merge annotated data with chromatine states with states with data
         crom_data <- merge(data, data.chrstates, by.x = "rs_number", by.y = "name" )

         fname <- paste0("ChrSates_Pla_data/List_CpGs_",
                         tools::file_path_sans_ext(basename(FilesToEnrich[i])),
                         "_annot_plac_chr_states.txt")
         dir.create("ChrSates_Pla_data", showWarnings = FALSE)
         write.table( crom_data, fname, quote=F, row.names=F, sep="\t")

         if("FDR" %in% colnames(data) | "Bonferroni" %in% colnames(data) & (BN==TRUE | !is.na(FDR) ))
         {

            if( tolower(testdata) =='fisher') {

               if( !is.na(FDR) ) {
                  ## --  HyperGeometric Test - States15_FP - BN,  BN_hyper and BN_hypo (Depletion and Enrichment)
                  States15FP_fdr <- getAllFisherTest(crom_data$bFDR , crom_data$States15_FP, outputdir = "ChrSates_15_Pla/Fisher_FDR", outputfile = FilesToEnrich[i])
                  States15FP_fdrhyper <- getAllFisherTest(FDR_Hyper, crom_data$States15_FP, outputdir = "ChrSates_15_Pla/Fisher_FDRHyper", outputfile = FilesToEnrich[i])
                  States15FP_fdrhypo <- getAllFisherTest(FDR_Hypo, crom_data$States15_FP, outputdir = "ChrSates_15_Pla/Fisher_FDRHypo", outputfile = FilesToEnrich[i])
               }

               if ( BN == TRUE) {
                  ## --  HyperGeometric Test - States15_FP - BN,  BN_hyper and BN_hypo (Depletion and Enrichment)
                  States15FP_bn <- getAllFisherTest(crom_data$Bonferroni, crom_data$States15_FP, outputdir = "ChrSates_15_Pla/Fisher_BN", outputfile = FilesToEnrich[i])
                  States15FP_bnhyper <- getAllFisherTest(BN_Hyper, crom_data$States15_FP, outputdir = "ChrSates_15_Pla/Fisher_BNHyper", outputfile = FilesToEnrich[i])
                  States15FP_bnhypo <- getAllFisherTest(BN_Hypo, crom_data$States15_FP, outputdir = "ChrSates_15_Pla/Fisher_BNHypo", outputfile = FilesToEnrich[i])
               }


            } else if ( tolower(testdata) =='hypergeometric') {

               if( !is.na(FDR) ) {
                  ## --  HyperGeometric Test - States15_FP - BN,  BN_hyper and BN_hypo (Depletion and Enrichment)
                  States15FP_fdr <- getAllHypergeometricTest(crom_data$FDR, crom_data$States15_FP, outputdir = "ChrSates_15_Pla/HyperG_FDR", outputfile = FilesToEnrich[i])
                  States15FP_fdrhyper <- getAllHypergeometricTest(FDR_Hyper, crom_data$States15_FP, outputdir = "ChrSates_15_Pla/HyperG_FDRHyper", outputfile = FilesToEnrich[i])
                  States15FP_fdrhypo <- getAllHypergeometricTest(FDR_Hypo, crom_data$States15_FP, outputdir = "ChrSates_15_Pla/HyperG_FDRHypo", outputfile = FilesToEnrich[i])

                  ## --  Resume in a table - HyperGeometric Test - States15_FP - FDR
                  resdata <- summary_States_FP_Table( crom_data$FDR, FDR_Hyper, FDR_Hypo, crom_data$States15_FP, outputdir = "ChrSates_15_Pla/Summary_HyperG_FDR", outputfile = FilesToEnrich[i], plot = TRUE )

                  ## --  Plot collapsed data HyperGeometric Test - States15_FP - FDR
                  plot_TestResults_Collapsed(list(fdr = hypergeo_States15FP_fdr, fdr_hypo = hypergeo_States15FP_fdrhypo, fdr_hyper = hypergeo_States15FP_fdrhyper),
                                             outputdir = "ChrSates_15_Pla", outputfile = FilesToEnrich[i])
               }

               if ( BN == TRUE) {
                  ## --  HyperGeometric Test - States15_FP - BN,  BN_hyper and BN_hypo (Depletion and Enrichment)
                  States15FP_bn <- getAllHypergeometricTest(crom_data$Bonferroni, crom_data$States15_FP, outputdir = "ChrSates_15_Pla/HyperG_BN", outputfile = FilesToEnrich[i])
                  States15FP_bnhyper <- getAllHypergeometricTest(BN_Hyper, crom_data$States15_FP, outputdir = "ChrSates_15_Pla/HyperG_BNHyper", outputfile = FilesToEnrich[i])
                  States15FP_bnhypo <- getAllHypergeometricTest(BN_Hypo, crom_data$States15_FP, outputdir = "ChrSates_15_Pla/HyperG_BNHypo", outputfile = FilesToEnrich[i])

                  ## --  Resume in a table - HyperGeometric Test - States15_FP - BN
                  resdata <- summary_States_FP_Table( crom_data$Bonferroni, BN_Hyper, BN_Hypo, crom_data$States15_FP, outputdir = "ChrSates_15_Pla/Summary_HyperG_BN", outputfile = FilesToEnrich[i], plot = TRUE )
               }
            }


            if ( !is.na(FDR) )  {
               ## --  Plot collapsed data HyperGeometric Test - States15_FP - FDR
               plot_TestResults_Collapsed(list(fdr = States15FP_fdr, fdr_hypo = States15FP_fdrhypo, fdr_hyper = States15FP_fdrhyper),
                                          outputdir = "ChrSates_15_Pla", outputfile = FilesToEnrich[i])
            }

            if ( BN == TRUE) {
               ## --  Plot collapsed data HyperGeometric Test - States15_FP - BN
               plot_TestResults_Collapsed(list(bn = States15FP_bn, bn_hypo = States15FP_bnhypo, bn_hyper = States15FP_bnhyper),
                                          outputdir = "ChrSates_15_Pla", outputfile = FilesToEnrich[i])
            }



            if(enrichFP18 == TRUE)
            {
               if( tolower(testdata) =='fisher') {

                  if( !is.na(FDR) ) {
                     ## --  HyperGeometric Test - States18_FP - BN,  BN_hyper and BN_hypo (Depletion and Enrichment)
                     States18FP_fdr <- getAllFisherTest(crom_data$bFDR , crom_data$States18_FP, outputdir = "ChrSates_18_Pla/Fisher_FDR", outputfile = FilesToEnrich[i])
                     States18FP_fdrhyper <- getAllFisherTest(FDR_Hyper, crom_data$States18_FP, outputdir = "ChrSates_18_Pla/Fisher_FDRHyper", outputfile = FilesToEnrich[i])
                     States18FP_fdrhypo <- getAllFisherTest(FDR_Hypo, crom_data$States18_FP, outputdir = "ChrSates_18_Pla/Fisher_FDRHypo", outputfile = FilesToEnrich[i])
                  }

                  if ( BN == TRUE) {
                     ## --  HyperGeometric Test - States18_FP - BN,  BN_hyper and BN_hypo (Depletion and Enrichment)
                     States18FP_bn <- getAllFisherTest(crom_data$Bonferroni, crom_data$States18_FP, outputdir = "ChrSates_18_Pla/Fisher_BN", outputfile = FilesToEnrich[i])
                     States18FP_bnhyper <- getAllFisherTest(BN_Hyper, crom_data$States18_FP, outputdir = "ChrSates_18_Pla/Fisher_BNHyper", outputfile = FilesToEnrich[i])
                     States18FP_bnhypo <- getAllFisherTest(BN_Hypo, crom_data$States18_FP, outputdir = "ChrSates_18_Pla/Fisher_BNHypo", outputfile = FilesToEnrich[i])
                  }

               } else if ( tolower(testdata) =='hypergeometric') {

                  if( !is.na(FDR) ) {
                     ## --  HyperGeometric Test - States18_FP - BN,  BN_hyper and BN_hypo (Depletion and Enrichment)
                     States18FP_fdr <- getAllHypergeometricTest(crom_data$FDR, crom_data$States18_FP, outputdir = "ChrSates_18_Pla/HyperG_FDR", outputfile = FilesToEnrich[i])
                     States18FP_fdrhyper <- getAllHypergeometricTest(FDR_Hyper, crom_data$States18_FP, outputdir = "ChrSates_18_Pla/HyperG_FDRHyper", outputfile = FilesToEnrich[i])
                     States18FP_fdrhypo <- getAllHypergeometricTest(FDR_Hypo, crom_data$States18_FP, outputdir = "ChrSates_18_Pla/HyperG_FDRHypo", outputfile = FilesToEnrich[i])

                     ## --  Resume in a table - HyperGeometric Test - States18_FP - BN
                     resdata <- summary_States_FP_Table( crom_data$FDR, FDR_Hyper, FDR_Hypo, crom_data$States18_FP, outputdir = "ChrSates_18_Pla/Summary_HyperG_FDR", outputfile = FilesToEnrich[i], plot = TRUE )
                  }

                  if ( BN == TRUE) {
                     ## --  HyperGeometric Test - States18_FP - BN,  BN_hyper and BN_hypo (Depletion and Enrichment)
                     States18FP_bn <- getAllHypergeometricTest(crom_data$Bonferroni, crom_data$States18_FP, outputdir = "ChrSates_18_Pla/HyperG_BN", outputfile = FilesToEnrich[i])
                     States18FP_bnhyper <- getAllHypergeometricTest(BN_Hyper, crom_data$States18_FP, outputdir = "ChrSates_18_Pla/HyperG_BNHyper", outputfile = FilesToEnrich[i])
                     States18FP_bnhypo <- getAllHypergeometricTest(BN_Hypo, crom_data$States18_FP, outputdir = "ChrSates_18_Pla/HyperG_BNHypo", outputfile = FilesToEnrich[i])

                     ## --  Resume in a table - HyperGeometric Test - States18_FP - BN
                     resdata <- summary_States_FP_Table( crom_data$Bonferroni, BN_Hyper, BN_Hypo, crom_data$States18_FP, outputdir = "ChrSates_18_Pla/Summary_HyperG_BN", outputfile = FilesToEnrich[i], plot = TRUE )
                  }
               }

               if ( !is.na(FDR) )  {
                  ## --  Plot collapsed data HyperGeometric Test - States15_FP - FDR
                  plot_TestResults_Collapsed(list(fdr = States18FP_fdr, fdr_hypo = States18FP_fdrhypo, fdr_hyper = States18FP_fdrhyper),
                                             outputdir = "ChrSates_18_Pla", outputfile = FilesToEnrich[i])
               }

               if ( BN == TRUE) {
                  ## --  Plot collapsed data HyperGeometric Test - States15_FP - BN
                  plot_TestResults_Collapsed(list(bn = States18FP_bn, bn_hypo = States18FP_bnhypo, bn_hyper = States18FP_bnhyper),
                                             outputdir = "ChrSates_18_Pla", outputfile = FilesToEnrich[i])
               }

            }
         }else {

            if( tolower(testdata) =='fisher') {
               States15FP <- getAllFisherTest(crom_data$signif , crom_data$States15_FP, outputdir = "ChrSates_15_Pla/Fisher_CpGlist", outputfile = FilesToEnrich[i])
            } else if ( tolower(testdata) =='hypergeometric') {
               States15FP <- getAllHypergeometricTest(crom_data$signif, crom_data$States15_FP, outputdir = "ChrSates_15_Pla/HyperG_CpGlist", outputfile = FilesToEnrich[i])
            }

            ## --  Plot collapsed data HyperGeometric Test - States15_FP - FDR
            #....CREC QUE NO ES NECESSARI.... # plot_OR(signif = States15FP, outputdir = "ChrSates_15_Pla", outputfile = FilesToEnrich[i])


            if(enrichFP18 == TRUE)
            {
               if( tolower(testdata) =='fisher') {
                  States18FP <- getAllFisherTest(crom_data$signif , crom_data$States18_FP, outputdir = "ChrSates_18_Pla/Fisher_CpGlist", outputfile = FilesToEnrich[i])
               } else if ( tolower(testdata) =='hypergeometric') {
                  States18FP <- getAllHypergeometricTest(crom_data$signif, crom_data$States18_FP, outputdir = "ChrSates_18_Pla/HyperG_CpGlist", outputfile = FilesToEnrich[i])
               }

               ## --  Plot collapsed data HyperGeometric Test - States15_FP - FDR
               #....CREC QUE NO ES NECESSARI.... # plot_TestResults_Collapsed(list(signif = States18FP), outputdir = "ChrSates_18_Pla", outputfile = FilesToEnrich[i])
            }

         }


         ## -- Partially Methylated Domains (PMDs) PLACENTA
         ## ------------------------------------------------

         # Create genomic ranges from PMD data
         PMD.GRange <- getEnrichGenomicRanges(PMD_placenta$Chr_PMD, PMD_placenta$Start_PMD, PMD_placenta$End_PMD)

         # Find overlaps between CpGs and PMD (find subject hits, query hits )
         overPMD <- findOverlapValues(data.GRange, PMD.GRange )

         #Create a data.frame with CpGs and PMDs information
         mdata <- cbind.data.frame(DataFrame(CpG = data.GRange$name[overPMD$qhits]), DataFrame(PMD = PMD.GRange$name[overPMD$shits]))

         # Merge with results from meta-analysis (A2)
         crom_data <- merge(crom_data, mdata, by.x="rs_number", by.y="CpG",all=T)
         # crom_data <- crom_data[order(crom_data$p.value),]

         # CpGs with PMD as NA
         PMD_NaN <- ifelse(is.na(crom_data$PMD),'IsNA','NotNA' )

         if("FDR" %in% colnames(data) | "Bonferroni" %in% colnames(data) & (BN==TRUE | !is.na(FDR) ))
         {

            if( tolower(testdata) =='fisher') {

               if( !is.na(FDR) ) {
                  ## --  HyperGeometric Test - PMD - FDR,  FDR_hyper and FDR_hypo  (Full data ) (Depletion and Enrichment)
                  PMD_fdr <- getAllFisherTest(crom_data$bFDR, PMD_NaN, outputdir = "PMD_Pla/Fisher_FDR", outputfile = FilesToEnrich[i])
                  PMD_fdrhyper <- getAllFisherTest(FDR_Hyper, PMD_NaN, outputdir = "PMD_Pla/Fisher_FDRHyper", outputfile = FilesToEnrich[i])
                  PMD_fdrhypo <- getAllFisherTest(FDR_Hypo, PMD_NaN, outputdir = "PMD_Pla/Fisher_FDRHypo", outputfile = FilesToEnrich[i])
               }
               if ( BN == TRUE) {
                  ## --  HyperGeometric Test - PMD - BN,  BN_hyper and BN_hypo  (Full data ) (Depletion and Enrichment)
                  PMD_bn <- getAllFisherTest(crom_data$Bonferroni, PMD_NaN, outputdir = "PMD_Pla/Fisher_BN", outputfile = FilesToEnrich[i])
                  PMD_bnhyper <- getAllFisherTest(BN_Hyper, PMD_NaN, outputdir = "PMD_Pla/Fisher_BNHyper", outputfile = FilesToEnrich[i])
                  PMD_bnhypo <- getAllFisherTest(BN_Hypo, PMD_NaN, outputdir = "PMD_Pla/Fisher_BNHypo", outputfile = FilesToEnrich[i])
               }

            } else if ( tolower(testdata) =='hypergeometric') {

               if( !is.na(FDR) ) {
                  ## --  HyperGeometric Test - PMD - FDR,  FDR_hyper and FDR_hypo  (Full data ) (Depletion and Enrichment)
                  PMD_fdr <- getAllHypergeometricTest(crom_data$bFDR, PMD_NaN, outputdir = "PMD_Pla/HyperG_FDR", outputfile = FilesToEnrich[i])
                  PMD_fdrhyper <- getAllHypergeometricTest(FDR_Hyper, PMD_NaN, outputdir = "PMD_Pla/HyperG_FDRHyper", outputfile = FilesToEnrich[i])
                  PMD_fdrhypo <- getAllHypergeometricTest(FDR_Hypo, PMD_NaN, outputdir = "PMD_Pla/HyperG_FDRHypo", outputfile = FilesToEnrich[i])
                  # Summary
                  resdata <- summary_HyperGeometrics_Table( crom_data$bFDR, FDR_Hyper, FDR_Hypo, PMD_NaN, outputdir = "PMD_Pla/Summary_HyperG_FDR", outputfile = FilesToEnrich[i], plot = TRUE )
               }

               if ( BN == TRUE) {
                  ## --  HyperGeometric Test - PMD - BN,  BN_hyper and BN_hypo  (Full data ) (Depletion and Enrichment)
                  PMD_bn <- getAllHypergeometricTest(crom_data$Bonferroni, PMD_NaN, outputdir = "PMD_Pla/HyperG_BN", outputfile = FilesToEnrich[i])
                  PMD_bnhyper <- getAllHypergeometricTest(BN_Hyper, PMD_NaN, outputdir = "PMD_Pla/HyperG_BNHyper", outputfile = FilesToEnrich[i])
                  PMD_bnhypo <- getAllHypergeometricTest(BN_Hypo, PMD_NaN, outputdir = "PMD_Pla/HyperG_BNHypo", outputfile = FilesToEnrich[i])
                  # Summary
                  resdata <- summary_HyperGeometrics_Table( crom_data$Bonferroni, BN_Hyper, BN_Hypo, PMD_NaN, outputdir = "PMD_Pla/Summary_HyperG_BN", outputfile = FilesToEnrich[i], plot = TRUE )
               }
            }

            if ( !is.na(FDR) )  {
               ## --  Plot collapsed data HyperGeometric Test - States15_FP - FDR
               plot_TestResults_Collapsed(list(fdr = PMD_fdr, fdr_hypo = PMD_fdrhypo, fdr_hyper = PMD_fdrhyper),
                                          outputdir = "PMD_Pla", outputfile = FilesToEnrich[i])
            }

            if ( BN == TRUE) {
               ## --  Plot collapsed data HyperGeometric Test - States15_FP - BN
               plot_TestResults_Collapsed(list(bn = PMD_bn, bn_hypo = PMD_bnhypo, bn_hyper = PMD_bnhyper),
                                          outputdir = "PMD_Pla", outputfile = FilesToEnrich[i])
            }
         } else {

            if( tolower(testdata) =='fisher') {
                  PMD <- getAllFisherTest(crom_data$signif, PMD_NaN, outputdir = "PMD_Pla/Fisher_CpGlist", outputfile = FilesToEnrich[i])
            } else if ( tolower(testdata) =='hypergeometric') {
                  PMD <- getAllHypergeometricTest(crom_data$signif, PMD_NaN, outputdir = "PMD_Pla/HyperG_CpGlist", outputfile = FilesToEnrich[i])
            }

         }



         ## -- Imprinting Regions PLACENTA
         ## ------------------------------------------------

         # Create genomic ranges from DMR data
         DMR.GRange <- getEnrichGenomicRanges(IR_Placenta$Chr_DMR, IR_Placenta$Start_DMR, IR_Placenta$End_DMR)

         # Find overlaps between CpGs and DMR (find subject hits, query hits )
         overDMR <- findOverlapValues(data.GRange, DMR.GRange )

         #Create a data.frame with CpGs and DMRs information
         mdata <- cbind.data.frame(DataFrame(CpG = data.GRange$name[overDMR$qhits]), DataFrame(DMR = DMR.GRange$name[overDMR$shits]))

         # Merge with results from meta-analysis (A2)
         crom_data <- merge(crom_data, mdata, by.x="rs_number", by.y="CpG",all=T)

         # CpGs with DMR as NA
         DMR_NaN <- ifelse(is.na(crom_data$DMR.y),'IsNA','NotNA' )

         if("FDR" %in% colnames(data) | "Bonferroni" %in% colnames(data) & (BN==TRUE | !is.na(FDR) ))
         {

            if( tolower(testdata) =='fisher') {

               if( !is.na(FDR) ) {
                  ## --  HyperGeometric Test - DMR - FDR,  FDR_hyper and FDR_hypo  (Full data ) (Depletion and Enrichment)
                  DMR_fdr <- getAllFisherTest(crom_data$bFDR, DMR_NaN, outputdir = "DMR_Pla/Fisher_FDR", outputfile = FilesToEnrich[i])
                  DMR_fdrhyper <- getAllFisherTest(FDR_Hyper, DMR_NaN, outputdir = "DMR_Pla/Fisher_FDRHyper", outputfile = FilesToEnrich[i])
                  DMR_fdrhypo <- getAllFisherTest(FDR_Hypo, DMR_NaN, outputdir = "DMR_Pla/Fisher_FDRHypo", outputfile = FilesToEnrich[i])
               }
               if ( BN == TRUE) {
                  ## --  HyperGeometric Test - DMR - BN,  BN_hyper and BN_hypo  (Full data ) (Depletion and Enrichment)
                  DMR_bn <- getAllFisherTest(crom_data$Bonferroni, DMR_NaN, outputdir = "DMR_Pla/Fisher_BN", outputfile = FilesToEnrich[i])
                  DMR_bnhyper <- getAllFisherTest(BN_Hyper, DMR_NaN, outputdir = "DMR_Pla/Fisher_BNHyper", outputfile = FilesToEnrich[i])
                  DMR_bnhypo <- getAllFisherTest(BN_Hypo, DMR_NaN, outputdir = "DMR_Pla/Fisher_BNHypo", outputfile = FilesToEnrich[i])
               }

            } else if ( tolower(testdata) =='hypergeometric') {

               if( !is.na(FDR) ) {
                  ## --  HyperGeometric Test - DMR - FDR,  FDR_hyper and FDR_hypo  (Full data ) (Depletion and Enrichment)
                  DMR_fdr <- getAllHypergeometricTest(crom_data$bFDR, DMR_NaN, outputdir = "DMR_Pla/HyperG_FDR", outputfile = FilesToEnrich[i])
                  DMR_fdrhyper <- getAllHypergeometricTest(FDR_Hyper, DMR_NaN, outputdir = "DMR_Pla/HyperG_FDRHyper", outputfile = FilesToEnrich[i])
                  DMR_fdrhypo <- getAllHypergeometricTest(FDR_Hypo, DMR_NaN, outputdir = "DMR_Pla/HyperG_FDRHypo", outputfile = FilesToEnrich[i])
                  # Summary
                  resdata <- summary_HyperGeometrics_Table( crom_data$bFDR, FDR_Hyper, FDR_Hypo, DMR_NaN, outputdir = "DMR_Pla/Summary_HyperG_FDR", outputfile = FilesToEnrich[i], plot = TRUE )
               }

               if ( BN == TRUE) {
                  ## --  HyperGeometric Test - DMR - BN,  BN_hyper and BN_hypo  (Full data ) (Depletion and Enrichment)
                  DMR_bn <- getAllHypergeometricTest(crom_data$Bonferroni, DMR_NaN, outputdir = "DMR_Pla/HyperG_BN", outputfile = FilesToEnrich[i])
                  DMR_bnhyper <- getAllHypergeometricTest(BN_Hyper, DMR_NaN, outputdir = "DMR_Pla/HyperG_BNHyper", outputfile = FilesToEnrich[i])
                  DMR_bnhypo <- getAllHypergeometricTest(BN_Hypo, DMR_NaN, outputdir = "DMR_Pla/HyperG_BNHypo", outputfile = FilesToEnrich[i])
                  # Summary
                  resdata <- summary_HyperGeometrics_Table( crom_data$Bonferroni, BN_Hyper, BN_Hypo, DMR_NaN, outputdir = "DMR_Pla/Summary_HyperG_BN", outputfile = FilesToEnrich[i], plot = TRUE )
               }
            }

            if ( !is.na(FDR) )  {
               ## --  Plot collapsed data HyperGeometric Test - States15_FP - FDR
               plot_TestResults_Collapsed(list(fdr = DMR_fdr, fdr_hypo = DMR_fdrhypo, fdr_hyper = DMR_fdrhyper),
                                          outputdir = "DMR_Pla", outputfile = FilesToEnrich[i])
            }

            if ( BN == TRUE) {
               ## --  Plot collapsed data HyperGeometric Test - States15_FP - BN
               plot_TestResults_Collapsed(list(bn = DMR_bn, bn_hypo = DMR_bnhypo, bn_hyper = DMR_bnhyper),
                                          outputdir = "DMR_Pla", outputfile = FilesToEnrich[i])
            }

         } else {

            if( tolower(testdata) =='fisher') {
                  DMR <- getAllFisherTest(crom_data$signif, DMR_NaN, outputdir = "DMR_Pla/Fisher_CpGlist", outputfile = FilesToEnrich[i])
            } else if ( tolower(testdata) =='hypergeometric') {
                  DMR <- getAllHypergeometricTest(crom_data$signif, DMR_NaN, outputdir = "DMR_Pla/HyperG_CpGlist", outputfile = FilesToEnrich[i])
            }
         }
      }

      # WRITE FINAL ENRICHMENT DATA
      write.table( crom_data, paste0( getwd(), "/",tools::file_path_sans_ext(basename(FilesToEnrich[i])),"_Enriched.csv" ) , quote=F, row.names=F, sep="\t")
   }

} else{
   print ("Error no data to enrich.")
}
