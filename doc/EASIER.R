## ----setup_knitr, include=FALSE-----------------------------------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE,
                      cache=TRUE)
library(knitr)
library(tools)

## ----textformat, include=FALSE, eval=TRUE-------------------------------------
colorize <- function(x, color) {
  if (knitr::is_latex_output()) {
    sprintf("\\textcolor{%s}{%s}", color, x)
  } else if (knitr::is_html_output()) {
    sprintf("<span style='color: %s;'>%s</span>", color, 
      x)
  } else x
}

## ----installDependences, eval=FALSE-------------------------------------------
#  # Install devtools
#  install.packages("devtools")
#  
#  # Install required packages
#  devtools::source_url("https://raw.githubusercontent.com/isglobal-brge/EASIER/HEAD/installer.R")
#  
#  # Install EASIER R package
#  devtools::install_github("isglobal-brge/EASIER@HEAD",  build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)

## ----loadMethyTools, eval=TRUE------------------------------------------------
library(EASIER)
library(readtext)

## ----qcworkflow, echo=FALSE, out.width='100%',  fig.align='center',  fig.cap="\\label{fig:qcworkflow}Quality control flowchart. This flowchart is used in the script under test folder to perform the quality control (QualityControl.R). The most important step in this workflow is the first step where we have to define the variables, if variables are well defined all the process is 'automatic' ", fig.pos='ht'----
include_graphics("imgs/workflows/QCWorkflow.png") 

## ----QC_varfiles--------------------------------------------------------------

files <- c('data/PROJ1_Cohort3_Model1_date_v2.txt',
           'data/PROJ1_Cohort3_Model2_date_v2.txt',
           'data/PROJ1_Cohort2_Plate_ModelA1_20170309.txt',
           'data/PROJ1_Cohort2_Plate_ModelA2_20170309.txt',
           'data/Cohort1_Model1_20170713.txt',
           'data/Cohort1_Model2_20170713.txt')

## ----QCVarres-----------------------------------------------------------------
# Result folder
results_folder <- 'QC_Results'

## ----QCVarPrefix--------------------------------------------------------------
# Prefixes for each file
prefixes <- c('PROJ1_Cohort3_A1', 'PROJ1_Cohort3_A2',
              'PROJ1_Cohort2_A1','PROJ1_Cohort2_A2', 
              'Cohort1_A1', 'Cohort1_A2')

## ----QCVarartype--------------------------------------------------------------
# Array type, used : EPIC or 450K
artype <- c('450K', 'EPIC', '450K', 'EPIC', '450K', 'EPIC')

## ----QCVarexclude-------------------------------------------------------------
# Parameters to exclude CpGs
exclude <-  c( 'MASK_sub30_copy', 
               'MASK_extBase', 
               'MASK_mapping', 
               'MASK_typeINextBaseSwitch', 
               'control_probes', 
               'Unrel_450_EPIC_blood', 
               'Sex')

## ----QCVarethnic--------------------------------------------------------------
ethnic <- c('EUR', 'EUR', 'SAS', 'SAS', 'EUR', 'EUR')

## ----QCVarN-------------------------------------------------------------------
N <- c(100, 100, 166, 166, 240, 240 )
n <- c(NA)

## ----QCCode1------------------------------------------------------------------

# Variable declaration to perform precision plot
medianSE <- numeric(length(files))
value_N <- numeric(length(files))
cohort_label <- character(length(files))

# Prepare output folder for results (create if not exists)
if(!dir.exists(file.path(getwd(), results_folder )))
   suppressWarnings(dir.create(file.path(getwd(), results_folder)))


# IMPORTANT FOR A REAL ANALYSIS :

# To show the execution flow we perform the analysis with only one data
# file. Normally, we have more than one data file to analyze, for that
# reason, we execute the code inside a loop and we follow the execution
# flow for each file defined in `files` 
# So we need to uncomment the for instruction and remove i <- 1 assignment.

# for ( i in 1:length(files) )
# {

   # we force i <- 1 to execute the analysis only for the first variable
   # for real data we have to remove this line
   i <- 1
   
   # Prepare output subfolder for cohort-model results (create if not exists)
   if(!dir.exists(file.path(getwd(), results_folder, prefixes[i] )))
      suppressWarnings(dir.create(file.path(getwd(), results_folder, prefixes[i])))

   # Creates an empty file to resume all data if an old file exist  removes
   # the file and creates a new one
   fResumeName <- paste0( file.path(getwd(), results_folder, prefixes[i]),"/",prefixes[i], "_descriptives.txt")
   if ( file.exists(fResumeName) ) {
      file.remove(fResumeName)
   }
   file.create(fResumeName)

   # Read data.
   cohort <- read.table(files[i], header = TRUE, as.is = TRUE)
   print(paste0("Cohort file : ",files[i]," - readed OK", sep = " "))

## ----QCCodeRead---------------------------------------------------------------
   
# Read data.
cohort <- read.table(files[i], header = TRUE, as.is = TRUE)
print(paste0("Cohort file : ",files[i]," - readed OK", sep = " "))

## ----QCCodeRemoveNA, eval=TRUE------------------------------------------------
# Remove rows with NA from data
cohort <- clean_NA_from_data(cohort)

## ----QCCodeDescriptives, eval=TRUE--------------------------------------------
# Descriptives - Before CpGs deletion #
descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, before = TRUE)


## ----QCCodeRemovedupli--------------------------------------------------------
# Remove duplicates
test_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )


## ----QCCodefilterLR-----------------------------------------------------------
# Remove cpGs with low representation
   # first, we test if colname_NforProbe and pcMissingSampes are defined
   if( !exists("colname_NforProbe") ) { colname_NforProbe <- NULL }
   if( !exists("pcMissingSamples") ) { pcMissingSamples <- NULL }

# Remove cpGs with low representation
cohort <- filterLowRepresentedCpGsinCohort(cohort, 
                                           colname_NforProbe, 
                                           pcMissingSamples, N[i], 
                                           fileresume = fResumeName )

## ----QCCodeexclCpGs-----------------------------------------------------------
# Exclude CpGs not meet conditions
   cohort <- exclude_CpGs(cohort, "probeID", 
                          exclude, 
                          ethnic = ethnic[i], 
                          filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), 
                          fileresume = fResumeName, artype = artype[i] )

## ----QCcodedesclast, eval=FALSE-----------------------------------------------
#  # Descriptives - After CpGs deletion #
#   descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"),
#                     fResumeName, before = FALSE )

## ----QCCodeAdjust-------------------------------------------------------------

# data before adjustment
head(cohort)

# Adjust data by Bonferroni and FDR
# Adjust data by Bonferroni and FDR
   cohort <- adjust_data(cohort, "P_VAL", 
                         bn=TRUE, 
                         fdr=TRUE, fResumeName, N[i] )

# data after adjustment
head(cohort)

## ----QCCodeWriteQData---------------------------------------------------------
# Write QC complete data to external file
write_QCData(cohort, paste0(results_folder, '/',prefixes[i], '/',prefixes[i]))

## ----QCcodedistrplotse, echo=TRUE, fig.align='center',  fig.cap="\\label{fig:QCcodedistrplotse}SE distribution plot", fig.pos='ht'----
   ## Visualization - Plots

   # Distribution plot
   plot_distribution(cohort$SE, 
                     main = paste('Standard Errors of', prefixes[i]), 
                     xlab = 'SE')

## ----QCcodedistrplotpval, echo=TRUE,  fig.align='center',  fig.cap="\\label{fig:QCcodedistrplotpval}p-value distribution plot", fig.pos='ht'----
   ## Visualization - Plots

   plot_distribution(cohort$P_VAL, 
                     main = paste('p-values of', prefixes[i]), 
                     xlab = 'p-value')

## ----QCcodeqqplot, echo=TRUE,   fig.align='center',  fig.cap="\\label{fig:QCcodeqqplot}QQ-plots", fig.pos='ht'----
   # QQ-plot.
   qqman::qq(cohort$P_VAL,
             main = sprintf('QQ-plot of %s (lambda = %f)', prefixes[i], 
                            lambda = get_lambda(cohort,"P_VAL")))

## ----QCCodeVolcanoplot, echo=TRUE, out.width='100%',  fig.align='center',  fig.cap="\\label{fig:QCCodeVolcanoplot}Volcano Plot", fig.pos='ht'----
   # Volcano plot.
   plot_volcano(cohort, "BETA", "P_VAL", 
                main=paste('Volcano plot of', prefixes[i]) )


## ----QCCodePrecisionP, echo=TRUE, eval=FALSE----------------------------------
#  plot_precisionp(precplot.data.n,
#                  paste(results_folder,  "precision_SE_N.png", sep='/'),
#                  main = "Subgroup Precision Plot -  1/median(SE) vs sqrt(n)")

## ----invchr17, echo=FALSE, out.width='100%',  fig.align='center',  fig.cap="\\label{fig:invchr17}Precision plot for 7 different datasets ", fig.pos='ht'----
include_graphics("imgs/QC/precision_SE_N.png") 

## ----QCCodeBetasBox, echo=TRUE, eval=FALSE------------------------------------
#  plot_betas_boxplot(betas.data,
#                     paste(results_folder, 'BETAS_BoxPlot.pdf', sep="/"))

## ----QCPlotBetasBox, echo=FALSE, out.width='85%',  fig.align='center',  fig.cap="\\label{fig:QCPlotBetasBox}Betas Boxplot plot for 10 different datasets ", fig.pos='ht'----
include_graphics("imgs/QC/BETAS_BoxPlot.png") 

## ----metaworkflow, echo=FALSE, out.width='100%', fig.align='center', fig.cap="\\label{fig:metaworkflow}Meta-analysis flowchart. The script MetaAnalysis.R, contains the code to perform the steps indicated in the flowchart.The only editing required by the researcher is defining the initial variables, which are specific to each study. The rest of the script is automatic and does not need any editing.", fig.pos='ht'----
include_graphics("imgs/workflows/MetaAnalysisWorkflow.png")

## ----meta_variables, eval=TRUE------------------------------------------------

## -- Variable definition for Meta-Analysis -- ##

# Array type, used : EPIC or 450K
artype <- '450K'

# Define data for each meta-analysis
metafiles <- list(
   'MetaA1' = c('Cohort1_A1','PROJ1_Cohort2_A1', 'PROJ1_Cohort3_A1' ),
   'MetaA2' = c('Cohort1_A2','PROJ1_Cohort2_A2', 'PROJ1_Cohort3_A2' ))

# Define maximum percent missing for each CpG
pcentMissing <- 0.8 # CpGs with precense lower than pcentMissing after EWAS
                    # meta-analysis will be deleted from the study.

# Paths with QCResults and path to store GWAMA results
results_folder <- 'QC_Results'
results_gwama <- '.'

# Venn diagrams ==> IMPORTANT : maximum 5 meta-analysis by venn diagram
venndiag_threshold <- 0.05

venn_diagrams <- list(
   c("MetaA1", "MetaA2" )
)

## -- End Variable definition for Meta-Analysis -- ##

## ----metaVarGwama, eval=FALSE-------------------------------------------------
#  
#  # GWAMA binary path  (GWAMA IsGlobal Server sw05 and sw06 installation)
#  gwama.dir <- paste0(Sys.getenv("HOME"),
#                      "/data/EWAS_metaanalysis/1_QC_results_cohorts/GWAMA/")
#  

## ----metaFolders, eval=TRUE---------------------------------------------------

## Create directory for GWAMA configuration files and GWAMA_Results 
## inside the defined results_gwama variable defined before.
if(!dir.exists(file.path(getwd(), 
                         paste(results_gwama, "GWAMA", sep="/") )))
   suppressWarnings(dir.create(file.path(getwd(), 
                                         paste(results_gwama, "GWAMA", sep="/"))))

## Create directory for GWAMA_Results
outputfolder <- paste0(results_gwama, "/GWAMA_Results")
if(!dir.exists(file.path(getwd(), outputfolder )))
   suppressWarnings(dir.create(file.path(getwd(), outputfolder)))

# We create a map file for GWAMA --> Used in Manhattan plots.
# We only need to indicate the array type
hapmapfile <- paste(results_gwama,"GWAMA", "hapmap.map" ,sep = "/")
generate_hapmap_file(artype, hapmapfile)



## ----metaFolders2, eval=FALSE-------------------------------------------------
#  
#  list.lowCpGs <- NULL
#  
#  # Create folder for a meta-analysis in GWAMA folder, here we
#  # store the GWAMA input files for each meta-analysis,
#  # We create one for complete meta-analysis
#  if(!dir.exists(file.path(getwd(),
#                           paste(results_gwama,"GWAMA", names(metafiles)[metf],
#                                 sep="/") )))
#     suppressWarnings(dir.create(file.path(getwd(),
#                                           paste(results_gwama,"GWAMA",
#                                                 names(metafiles)[metf],
#                                                 sep="/"))))
#  # We create another for meta-analysis without filtered CpGs with low
#  # percentage (sufix _Filtr)
#  if(!dir.exists(file.path(getwd(),
#                           paste0(results_gwama,"/GWAMA/",
#                                  names(metafiles)[metf],
#                                  "_Filtr") )))
#     suppressWarnings(dir.create(file.path(getwd(),
#                                           paste0(results_gwama, "/GWAMA/",
#                                                  names(metafiles)[metf],
#                                                  "_Filtr"))))
#  
#  # GWAMA File name base
#  inputfolder <- paste0(results_gwama,"/GWAMA/",  names(metafiles)[metf])
#  
#  modelfiles <- unlist(metafiles[metf])
#  
#  # Execution with all CpGs and without filtered CpGs
#  runs <- c('Normal', 'lowcpgs')
#  lowCpGs = FALSE;
#  outputfiles <- list()
#  
#  outputgwama <- paste(outputfolder,names(metafiles)[metf],sep = '/')
#  

## ----metaAnalys, eval=FALSE---------------------------------------------------
#  
#   if(runs[j]=='lowcpgs') {
#     lowCpGs = TRUE
#     # Get low presence CpGs in order to exclude this from the new meta-analysis
#     list.lowCpGs <- get_low_presence_CpGs(outputfiles[[j-1]], pcentMissing)
#     inputfolder <- paste0(results_gwama,"/GWAMA/",  names(metafiles)[metf], "_Filtr")
#     outputgwama <- paste0(outputgwama,"_Filtr")
#   }
#  
#   # Create a GWAMA files for each file in meta-analysis and one file with
#   # gwama meta-analysis configuration
#   for ( i in 1:length(modelfiles) )
#     create_GWAMA_files(results_folder,  modelfiles[i],
#                        inputfolder, N[i], list.lowCpGs )

## ----metaGWANA, eval=FALSE----------------------------------------------------
#   # Execute GWAMA meta-analysis and manhattan-plot, QQ-plot and a file
#   # with gwama results.
#   outputfiles[[runs[j]]] <- run_GWAMA_MetaAnalysis(inputfolder,
#                                                   outputgwama,
#                                                   names(metafiles)[metf],
#                                                   gwama.dir)

## ----EnrichCommon_2_1, eval=TRUE----------------------------------------------
if(dim(data)[1] <= 1 | dim(data)[2] <= 1) {
   data <- read.table(FilesToEnrich[i], dec = ".") # Avoid header
   data <- as.vector(t(data))
   head(data)
   data <- get_annotattions(data, artype, FilesToEnrich[i], outputfolder )
   
   head(data)
   
   allCpGs <- TRUE
   data$chromosome <- substr(data$chr,4,length(data$chr))
   data$rs_number <- data$CpGs
}


## ----EnrichUniqueGene, eval=FALSE---------------------------------------------
#        # get unique genes from data
#  geneUniv <- lapply( lapply(miss_enrich[grepl("signif", names(miss_enrich))],
#                             function(cpgs) {
#                                data[which(as.character(data$CpGs) %in% cpgs),]$UCSC_RefGene_Name
#                                }),
#                      getUniqueGenes)
#  
#  geneUniv

## ----EnrichCommon_2, eval=FALSE, cache=TRUE-----------------------------------
#  ## -- Functional Enrichmnet
#  ## ------------------------
#  
#  # Enrichment with missMethyl - GO and KEGG --> Writes results to outputfolder
#  miss_enrich <- missMethyl_enrichment(data, outputfolder, FilesToEnrich[i],
#                                       artype, BN, FDR, pvalue, allCpGs, plots = TRUE )
#  
#  head(miss_enrich$GO)
#  head(miss_enrich$KEGG)

## ----EnrichCommon_3, eval=FALSE, cache=TRUE-----------------------------------
#  ## -- Molecular Enrichmnet
#  ## -----------------------
#  
#  # Molecular Signatures Database enrichment
#  msd_enrich <- MSigDB_enrichment(data, outputfolder, FilesToEnrich[i], artype, BN, FDR, pvalue, allCpGs)
#  
#  head(msd_enrich$MSigDB)
#  

## ----ConsensusPathDB, eval=FALSE, cache=TRUE----------------------------------
#  
#  ## -- Online Tools
#  
#  # Enrichment with ConsensusPathDB
#  #     - Consensus path http://consensuspathdb.org/
#  #     (gene-set analysis – over-representation analysis)
#  
#  # Available FSet types :
#  # 1 P     manually curated pathways from pathway databases
#  # 2 N     interaction network neighborhood-based functional sets
#  # 3 G2    Gene Ontology-based sets, GO level 2
#  # 4 G3    Gene Ontology-based sets, GO level 3
#  # 5 G4    Gene Ontology-based sets, GO level 4
#  # 6 G5    Gene Ontology-based sets, GO level 5
#  # 7 C     protein complex-based sets
#  
#  acFSet <- c('C', 'P', 'G2', 'G3')
#  acType <- 'entrez-gene'
#  
#  # Get Enrichment
#  CPDB_enrich <- lapply(names(geneUniv), function( data, accFSet, genes ) {
#     print(data)
#     lapply(accFSet,
#            get_consensusPdb_OverRepresentation,
#            entityType='genes',
#            accNumbers=na.omit(as.character(eval(parse(text = paste0("genes$",data))))),
#            accType=acType,
#            outputdir = "ConsensusPathDB",
#            outputfile = gsub(".", "_", data, fixed=TRUE) )},
#     accFSet = acFSet, genes = geneUniv)
#  
#  names(CPDB_enrich) <- names(geneUniv)
#  

## ----EnrichComStatsPrep, eval=FALSE-------------------------------------------
#  if("FDR" %in% colnames(data) & "Bonferroni" %in% colnames(data))
#  {
#  
#     ## -- Prepare data
#     ## ---------------
#  
#     # Add column bFDR to data for that CpGs that accomplish with FDR
#      # Classify fdr into "yes" and no taking into account FDR significance level
#     data$bFDR <- getBinaryClassificationYesNo(data$FDR, "<", FDR)
#  
#     # Classify by Hyper and Hypo methylated
#     data$meth_state <- getHyperHypo(data$beta) # Classify methylation into Hyper and Hypo
#  
#     # CpGs FDR and Hyper and Hypo respectively
#     FDR_Hyper <- ifelse(data$bFDR == 'yes' &
#                            data$meth_state=='Hyper', "yes", "no")
#     FDR_Hypo <- ifelse(data$bFDR == 'yes' &
#                           data$meth_state=='Hypo', "yes", "no")
#  
#     # CpGs Bonferroni and Hyper and Hypo respectively
#     BN_Hyper <- ifelse(data$Bonferroni == 'yes' &
#                           data$meth_state=='Hyper', "yes", "no")
#     BN_Hypo <- ifelse(data$Bonferroni == 'yes' &
#                          data$meth_state=='Hypo', "yes", "no")

## ----EnrichComStats_GeneP, eval=FALSE-----------------------------------------
#  ## --  CpG Gene position
#  ## ---------------------
#  
#  # Get descriptives
#  get_descriptives_GenePosition(data$UCSC_RefGene_Group,
#                                data$Bonferroni,
#                                "Bonferroni",
#                                outputdir = "GenePosition/Fisher_BN_Desc",
#                                outputfile = FilesToEnrich[i])
#  get_descriptives_GenePosition(data$UCSC_RefGene_Group, d
#                                ata$bFDR , "FDR",
#                                outputdir = "GenePosition/Fisher_FDR_Desc",
#                                outputfile = FilesToEnrich[i])
#  
#  
#  if( tolower(testdata) =='fisher') {
#     ## --  Fisher Test - Gene position - FDR, FDR_hyper and FDR_hypo
#     GenePosition_fdr <- getAllFisherTest(data$bFDR,
#                                    data$UCSC_RefGene_Group,
#                                    outputdir = "GenePosition/Fisher_FDR",
#                                    outputfile = FilesToEnrich[i],
#                                    plots = TRUE )
#     GenePosition_fdr_hyper <- getAllFisherTest(FDR_Hyper,
#                                    data$UCSC_RefGene_Group,
#                                    outputdir = "GenePosition/Fisher_FDRHyper",
#                                    outputfile = FilesToEnrich[i],
#                                    plots = TRUE )
#     GenePosition_fdr_hypo <- getAllFisherTest(FDR_Hypo,
#                                      data$UCSC_RefGene_Group,
#                                      outputdir = "GenePosition/Fisher_FDRHypo",
#                                      outputfile = FilesToEnrich[i], plots = TRUE )
#  }
#  else if ( tolower(testdata) =='hypergeometric') {
#     ## --  HyperGeometric Test - Island relative position -
#     ## FDR, FDR_hyper and FDR_hypo (for Depletion and Enrichment)
#     GenePosition_fdr <- getAllHypergeometricTest(data$bFDR,
#                                      data$UCSC_RefGene_Group,
#                                      outputdir = "GenePosition/HyperG_FDR",
#                                      outputfile = FilesToEnrich[i])
#     GenePosition_fdr_hyper <- getAllHypergeometricTest(FDR_Hyper,
#                                      data$UCSC_RefGene_Group,
#                                      outputdir = "GenePosition/HyperG_FDRHyper",
#                                      outputfile = FilesToEnrich[i])
#     GenePosition_fdr_hypo <- getAllHypergeometricTest(FDR_Hypo,
#                                     data$UCSC_RefGene_Group,
#                                     outputdir = "GenePosition/HyperG_FDRHypo",
#                                     outputfile = FilesToEnrich[i])
#  }

## ----EnrichComStatsGenePplots, eval=FALSE-------------------------------------
#           plot_TestResults_Collapsed(list(fdr = GenePosition_fdr,
#                                           fdr_hypo = GenePosition_fdr_hypo,
#                                           fdr_hyper = GenePosition_fdr_hyper),
#                                      outputdir = "GenePosition",
#                                      outputfile = FilesToEnrich[i], main = )

## ----enriGenePos, echo=FALSE, out.width='110%',  fig.align='center',  fig.cap="\\label{fig:genepos}Gene position with Fisher test for Hyper and Hypo methylated CpGs", fig.pos='ht'----
include_graphics("imgs/enrich/plot_genepos.png") 

## ----EnrichComStatsRelIsland, eval=FALSE--------------------------------------
#  
#  ## --  CpG Island relative position
#  ## --------------------------------
#  
#  # Get descriptives
#  get_descriptives_RelativetoIsland(data$Relation_to_Island,
#                              data$Bonferroni,
#                              "Bonferroni",
#                              outputdir = "RelativeToIsland/Fisher_BN_RelativeToIsland",
#                              outputfile = FilesToEnrich[i])
#  get_descriptives_RelativetoIsland(data$Relation_to_Island,
#                              data$bFDR ,
#                              "FDR",
#                              outputdir = "RelativeToIsland/Fisher_FDR_RelativeToIsland",
#                              outputfile = FilesToEnrich[i])
#  
#  
#  if( tolower(testdata) =='fisher') {
#     ## --  Fisher Test - Position Relative to Island - FDR, FDR_hyper and FDR_hypo
#     relative_island_fdr <- getAllFisherTest(data$bFDR,
#                                    data$Relation_to_Island,
#                                    outputdir = "RelativeToIsland/Fisher_FDR",
#                                    outputfile = FilesToEnrich[i], plots = TRUE )
#     relative_island_fdr_hyper <- getAllFisherTest(FDR_Hyper,
#                                    data$Relation_to_Island,
#                                    outputdir = "RelativeToIsland/Fisher_FDRHyper",
#                                    outputfile = FilesToEnrich[i], plots = TRUE )
#     relative_island_fdr_hypo <- getAllFisherTest(FDR_Hypo,
#                                         data$Relation_to_Island,
#                                         outputdir = "RelativeToIsland/Fisher_FDRHypo",
#                                         outputfile = FilesToEnrich[i], plots = TRUE )
#  }
#  
#  

## ----enrichRelIsland, echo=FALSE, out.width='110%',  fig.align='center',  fig.cap="\\label{fig:genepos}Gene position with Fisher test for Hyper and Hypo methylated CpGs", fig.pos='ht'----
include_graphics("imgs/enrich/plot_relatisland.png") 

## ----enrichworkflowblood, echo=FALSE, out.width='110%',  fig.align='center',  fig.cap="\\label{fig:enrichworkflowblood}Enrichment flowchart. Detailed Blood enrichment", fig.pos='ht'----
include_graphics("imgs/enrich/blood.png") 

## ----EnrichBlood, eval=FALSE--------------------------------------------------
#  
#  ## --  ROADMAP  -  Metilation in Cromatine States - BLOOD
#  ## -------------------------------------------------------
#  ##  Analysis of methylation changes in the different chromatin
#  ##  states (CpGs are diff meth in some states and others don't)
#  
#  # Prepare data
#  # Adds chromatine state columns
#  crom_data <- addCrom15Columns(data, "CpGId")
#  
#  if("FDR" %in% colnames(data) & "Bonferroni" %in% colnames(data))
#  {
#  
#     # Columns with chromatin status information :
#     ChrStatCols <- c("TssA","TssAFlnk","TxFlnk","TxWk","Tx","EnhG",
#                      "Enh","ZNF.Rpts","Het","TssBiv","BivFlnk",
#                      "EnhBiv","ReprPC","ReprPCWk","Quies")
#  
#     if( !is.na(FDR) ) {
#        chrom_states_fdr <- getAllChromStateOR( crom_data$bFDR,
#                                    crom_data[,ChrStatCols],
#                                    outputdir = "CromStates/OR_FDR",
#                                    outputfile = FilesToEnrich[i],
#                                    plots = TRUE )
#        chrom_states_fdr_hyper <- getAllChromStateOR( FDR_Hyper,
#                                    crom_data[,ChrStatCols],
#                                    outputdir = "CromStates/OR_FDRHyper",
#                                    outputfile = FilesToEnrich[i],
#                                    plots = TRUE )
#        chrom_states_fdr_hypo <- getAllChromStateOR( FDR_Hypo,
#                                   crom_data[,ChrStatCols],
#                                   outputdir = "CromStates/OR_FDRHypo",
#                                   outputfile = FilesToEnrich[i],
#                                   plots = TRUE )
#     }
#     if ( BN == TRUE) {
#        chrom_states_bn <- getAllChromStateOR( crom_data$Bonferroni,
#                                   crom_data[,ChrStatCols],
#                                   outputdir = "CromStates/OR_BN",
#                                   outputfile = FilesToEnrich[i],
#                                   plots = TRUE )
#        chrom_states_bn_hyper <- getAllChromStateOR( BN_Hyper,
#                                   crom_data[,ChrStatCols],
#                                   outputdir = "CromStates/OR_BNHyper",
#                                   outputfile = FilesToEnrich[i],
#                                   plots = TRUE )
#        chrom_states_bn_hypo <- getAllChromStateOR( BN_Hypo,
#                                  crom_data[,ChrStatCols],
#                                  outputdir = "CromStates/OR_BNHypo",
#                                  outputfile = FilesToEnrich[i],
#                                        plots = TRUE )
#     }
#  }
#  

## ----enrichworkflowplacenta, echo=FALSE, out.width='110%',  fig.align='center',  fig.cap="\\label{fig:enrichworkflowplacenta}Enrichment flowchart. Detailed Placenta enrichment", fig.pos='ht'----
include_graphics("imgs/enrich/placenta.png") 

## ----EnrichPlacentachromstates, eval=FALSE------------------------------------
#  
#  ## -- ROADMAP  -  Regulatory feature enrichment analysis - PLACENTA
#  ## -----------------------------------------------------------------
#  
#  # Convert to Genomic Ranges
#  data.GRange <- GRanges(
#     seqnames = Rle(data$chr),
#     ranges=IRanges(data$pos, end=data$pos),
#     name=data$CpGs,
#     chr=data$chromosome,
#     pos=data$pos
#  )
#  names(data.GRange) <- data.GRange$name
#  
#  # Find overlaps between CpGs and Fetal Placenta (States 15 and 18)
#  over15 <- findOverlapValues(data.GRange, FP_15_E091 )
#  
#  if (enrichFP18 == TRUE){
#     over18 <- findOverlapValues(data.GRange, FP_18_E091 )
#     # Add states 15 and 18 to data.GRange file
#     # and write to a file : CpGs, state15 and state18
#     data.chrstates <- c(mcols(over15$ranges), over15$values, over18$values)
#     colnames(data.chrstates)[grep("States",colnames(data.chrstates))] <-
#        c("States15_FP", "States18_FP")
#  } else {
#     # Add states 15 to data.GRange file and write to a file : CpGs, state15
#     data.chrstates <- c(mcols(over15$ranges), over15$values)
#     colnames(data.chrstates)[grep("States",colnames(data.chrstates))] <-
#        c("States15_FP")
#  }
#  
#  # Merge annotated data with chromatine states with states with data
#  crom_data <- merge(data, data.chrstates, by.x = "CpGs", by.y = "name" )
#  
#  fname <- paste0("ChrSates_Pla_data/List_CpGs_",
#                  tools::file_path_sans_ext(basename(FilesToEnrich[i])),
#                  "_annot_plac_chr_states.txt")
#  dir.create("ChrSates_Pla_data", showWarnings = FALSE)
#  write.table( crom_data, fname, quote=F, row.names=F, sep="\t")
#  
#  ## --  Fisher Test - States15_FP - BN,  BN_hyper and BN_hypo
#  ## (Depletion and Enrichment)
#  States15FP_bn <- getAllFisherTest(crom_data$Bonferroni,
#                                 crom_data$States15_FP,
#                                 outputdir = "ChrSates_15_Pla/Fisher_BN",
#                                 outputfile = FilesToEnrich[i])
#  States15FP_bnhyper <- getAllFisherTest(BN_Hyper,
#                                   crom_data$States15_FP,
#                                   outputdir = "ChrSates_15_Pla/Fisher_BNHyper",
#                                   outputfile = FilesToEnrich[i])
#  States15FP_bnhypo <- getAllFisherTest(BN_Hypo,
#                                  crom_data$States15_FP,
#                                  outputdir = "ChrSates_15_Pla/Fisher_BNHypo",
#                                  outputfile = FilesToEnrich[i])
#  
#  
#  ## --  Plot collapsed data HyperGeometric Test - States15_FP - BN
#  plot_TestResults_Collapsed(list(bn = States15FP_bn,
#                                  bn_hypo = States15FP_bnhypo,
#                                  bn_hyper = States15FP_bnhyper),
#                             outputdir = "ChrSates_15_Pla",
#                             outputfile = FilesToEnrich[i])
#  
#  

## ----chr15pla1, echo=FALSE, out.width='110%',  fig.align='center',  fig.cap="\\label{fig:chr15pla1}Chromatine states 15 for placenta - Fisher test for Hyper and Hypo methylated CpGs", fig.pos='ht'----
include_graphics("imgs/enrich/chr15pla.png") 

## ----EnrichPlacentaPMD, eval=FALSE--------------------------------------------
#  
#  ## -- Partially Methylated Domains (PMDs) PLACENTA
#  ## ------------------------------------------------
#  
#  # Create genomic ranges from PMD data
#  PMD.GRange <- getEnrichGenomicRanges(PMD_placenta$Chr_PMD,
#                                       PMD_placenta$Start_PMD,
#                                       PMD_placenta$End_PMD)
#  
#  # Find overlaps between CpGs and PMD (find subject hits, query hits )
#  overPMD <- findOverlapValues(data.GRange, PMD.GRange )
#  
#  #Create a data.frame with CpGs and PMDs information
#  mdata <- as.data.frame(cbind(DataFrame(CpG = data.GRange$name[overPMD$qhits]),
#                               DataFrame(PMD = PMD.GRange$name[overPMD$shits])))
#  
#  # Merge with results from meta-analysis (A2)
#  crom_data <- merge(crom_data, mdata, by.x="CpGs", by.y="CpG",all=T)
#  # crom_data <- crom_data[order(crom_data$p.value),]
#  
#  # CpGs with PMD as NA
#  PMD_NaN <- ifelse(is.na(crom_data$PMD),'IsNA','NotNA' )
#  
#  
#  ## --  Fisher Test - PMD - BN,  BN_hyper and BN_hypo
#  ## (Full data ) (Depletion and Enrichment)
#  PMD_bn <- getAllFisherTest(crom_data$Bonferroni,
#                             PMD_NaN,
#                             outputdir = "PMD_Pla/Fisher_BN",
#                             outputfile = FilesToEnrich[i])
#  PMD_bnhyper <- getAllFisherTest(BN_Hyper,
#                                  PMD_NaN,
#                                  outputdir = "PMD_Pla/Fisher_BNHyper",
#                                  outputfile = FilesToEnrich[i])
#  PMD_bnhypo <- getAllFisherTest(BN_Hypo,
#                                 PMD_NaN,
#                                 outputdir = "PMD_Pla/Fisher_BNHypo",
#                                 outputfile = FilesToEnrich[i])
#  
#   ## --  Plot collapsed data HyperGeometric Test - States15_FP - BN
#  plot_TestResults_Collapsed(list(bn = PMD_bn,
#                                  bn_hypo = PMD_bnhypo,
#                                  bn_hyper = PMD_bnhyper),
#                             outputdir = "PMD_Pla",
#                             outputfile = FilesToEnrich[i])
#  

## ----pmd15pla, echo=FALSE, out.width='110%',  fig.align='center',  fig.cap="\\label{fig:pmd15pla}Partial Metilated Domains for placenta - Fisher test for Hyper and Hypo methylated CpGs", fig.pos='ht'----
include_graphics("imgs/enrich/pmd_pla.png") 

## ----EnrichPlacentaIR, eval=FALSE---------------------------------------------
#  ## -- Imprinting Regions PLACENTA
#  ## -------------------------------
#  
#  # Create genomic ranges from DMR data
#  DMR.GRange <- getEnrichGenomicRanges(IR_Placenta$Chr_DMR,
#                                       IR_Placenta$Start_DMR,
#                                       IR_Placenta$End_DMR)
#  
#  # Find overlaps between CpGs and DMR (find subject hits, query hits )
#  overDMR <- findOverlapValues(data.GRange, DMR.GRange )
#  
#  #Create a data.frame with CpGs and DMRs information
#  mdata <- as.data.frame(cbind(DataFrame(CpG = data.GRange$name[overDMR$qhits]),
#                               DataFrame(DMR = DMR.GRange$name[overDMR$shits])))
#  
#  # Merge with results from meta-analysis (A2)
#  crom_data <- merge(crom_data, mdata, by.x="CpGs", by.y="CpG",all=T)
#  
#  # CpGs with DMR as NA
#  DMR_NaN <- ifelse(is.na(crom_data$DMR.y),'IsNA','NotNA' )
#  
#  ## --  Fisher Test - DMR - BN,  BN_hyper and BN_hypo
#  ## (Full data ) (Depletion and Enrichment)
#  DMR_bn <- getAllFisherTest(crom_data$Bonferroni,
#                             DMR_NaN,
#                             outputdir = "DMR_Pla/Fisher_BN",
#                             outputfile = FilesToEnrich[i])
#  DMR_bnhyper <- getAllFisherTest(BN_Hyper,
#                                  DMR_NaN,
#                                  outputdir = "DMR_Pla/Fisher_BNHyper",
#                                  outputfile = FilesToEnrich[i])
#  DMR_bnhypo <- getAllFisherTest(BN_Hypo,
#                                 DMR_NaN,
#                                 outputdir = "DMR_Pla/Fisher_BNHypo",
#                                 outputfile = FilesToEnrich[i])
#  
#   ## --  Plot collapsed data HyperGeometric Test - States15_FP - BN
#  plot_TestResults_Collapsed(list(bn = DMR_bn,
#                                  bn_hypo = DMR_bnhypo,
#                                  bn_hyper = DMR_bnhyper),
#                             outputdir = "DMR_Pla", outputfile = FilesToEnrich[i])
#  

## ----chr15pla, echo=FALSE, out.width='110%',  fig.align='center',  fig.cap="\\label{fig:chr15pla}Imprinted Regions for placenta - Fisher test for Hyper and Hypo methylated CpGs", fig.pos='ht'----
include_graphics("imgs/enrich/iregions.png") 

