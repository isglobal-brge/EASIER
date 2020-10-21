if (!require(rasterpdf, quietly = TRUE)) install.packages('rasterpdf', repos = 'https://cran.rediris.es/' )
if (!require(meta, quietly = TRUE)) install.packages('meta', repos = 'https://cran.rediris.es/' )
if (!require(tibble, quietly = TRUE)) install.packages('tibble')
if (!require(dplyr, quietly = TRUE)) install.packages('dplyr')
if (!require(tidyverse, quietly = TRUE)) install.packages::install( "tidyverse" )
if (!require(stringr, quietly = TRUE)) install.packages('stringr')
if (!require(meta, quietly = TRUE)) install.packages('meta') # Forest Plot
if (!require(missMethyl, quietly = TRUE)) BiocManager::install( "missMethyl" )
if (!require(org.Hs.eg.db, quietly = TRUE)) BiocManager::install( "org.Hs.eg.db" )
if (!require(GenomicRanges, quietly = TRUE)) BiocManager::install( "GenomicRanges" )
if (!require(rtracklayer, quietly = TRUE)) BiocManager::install( "rtracklayer" )


# Plots
if (!require(ggplot2, quietly = TRUE)) install.packages('ggplot2')
if (!require(VennDiagram, quietly = TRUE)) install.packages('VennDiagram')
if (!require(RColorBrewer, quietly = TRUE)) install.packages('RColorBrewer')
if (!require(reshape, quietly = TRUE)) install.packages('reshape')
if (!require(ggsignif, quietly = TRUE)) install.packages('ggsignif')


if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")

#BiocManager::install( c("IlluminaHumanMethylation450kanno.ilmn12.hg19",
#                        "IlluminaHumanMethylation450kanno.ilmn12.hg19",
#                        "missMethyl",
#                        "org.Hs.eg.db",
#                        "GenomicRanges") )


library(rasterpdf)
library(meta)
library(tibble)
library(dplyr)
library(stringr)
library(meta)
library(missMethyl)
library(org.Hs.eg.db)
library(GenomicRanges)
library(rtracklayer)

# Plots
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
library(reshape)
library(ggsignif)






library(methyTools)

#..# setwd("~/Library/Mobile Documents/com~apple~CloudDocs/PROJECTES/Treballant/methyTools")
#..# setwd("/Users/mailos/Library/Mobile Documents/com~apple~CloudDocs/ISGlobal/Projectes/EWAS/2_Llibreria/proves")
setwd("/Users/mailos/tmp/proves")


# Study name
study <- 'TEST'

files <- c('data/PACE_AQUA_Model1_date_v2.txt',
           'data/PACE_AQUA_Model2_date_v2.txt',
           'data/PACE_INMA_Plate_ModelA1_20170309.txt',
           'data/PACE_INMA_Plate_ModelA2_20170309.txt',
           'data/PACE_INMA_Plate_ModelB1_20170320.txt',
           'data/PACE_INMA_Plate_ModelB2_20170320.txt',
           'data/PACE_INMA_Plate_ModelC1_20170818.txt',
           'data/PACE_INMA_Plate_ModelC2_20170818.txt',
           'data/RICHS_Model1_20170713.txt',
           'data/RICHS_Model2_20170713.txt')

# Result folder
results_folder <- 'QC_Results'

# Prefixes for each file
prefixes <- c('PACE_AQUA_A1', 'PACE_AQUA_A2',
              'PACE_IMMA_A1','PACE_IMMA_A2', 'PACE_IMMA_B1', 'PACE_IMMA_B2', 'PACE_IMMA_C1', 'PACE_IMMA_C2',
              'RICHS_A1', 'RICHS_A2')

# Array type, used : EPIC or 450K
artype <- '450K'

# Parameters to exclude CpGs
exclude <- c( 'MASK_sub35_copy', 'MASK_typeINextBaseSwitch', 'noncpg_probes', 'control_probes', 'Unreliable_450_EPIC', 'Sex')

ethnic <- 'EUR'

N <- c(100, 100, 166, 166, 166, 166, 166, 166, 240, 240 )
n <- c(NA)


# Venn diagrams
venn_diagrams <- list(
   c("PACE_AQUA_A1", "PACE_IMMA_A1", "PACE_IMMA_B1", "PACE_IMMA_C1", "RICHS_A1" ),
   c("PACE_AQUA_A2", "PACE_IMMA_A2", "PACE_IMMA_B2", "PACE_IMMA_C2", "RICHS_A2" )
)



## ###################### ##
##  QC - Quality Control  ##
## ###################### ##

# Variable declaration to perform precision plot
medianSE <- numeric(length(files))
value_N <- numeric(length(files))
cohort_label <- character(length(files))

# Prepare output folder for results (create if not exists)
if(!dir.exists(file.path(getwd(), results_folder )))
   suppressWarnings(dir.create(file.path(getwd(), results_folder)))

## Remove duplicates, Exclude CpGs and adjust data (BN and FDR)
for ( i in 1:length(files) )
{

   # Read data.
   cohort <- read.table(files[i], header = TRUE, as.is = TRUE)
   print(paste0("Cohort file : ",files[i]," - readed OK", sep = " "))

   # Remove rows with NA from data
   cohort <- clean_NA_from_data(cohort)

   # Descriptives - Before CpGs deletion #
   descriptives_CpGs(cohort, seq(2,4), paste0(results_folder,'/',prefixes[i],'_descriptives_init.txt') )

   # Remove duplicates
   cohort <- remove_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i],'_descriptives_duplic.txt'), paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )

   # Exclude CpGs not meet conditions
   cohort <- exclude_CpGs(cohort, "probeID", exclude, filename = paste0(results_folder,'/',prefixes[i],'_excluded.txt') )

   # Descriptives - After CpGs deletion #
   descriptives_CpGs(cohort, seq(2,4), paste0(results_folder,'/',prefixes[i],'_descriptives_last.txt') )

   # Adjust data by Bonferroni and FDR
   cohort <- adjust_data(cohort, "P_VAL", bn=TRUE, fdr=TRUE, filename =  paste0(results_folder,'/',prefixes[i],'_ResumeSignificatives.txt')  )

   # Write QC complete data to external file
   write_QCData(cohort, paste0(results_folder,'/',prefixes[i]))


   ## Visualization - Plots
   rasterpdf::raster_pdf(paste0(results_folder,'/',prefixes[i], '_QCplots.pdf'), res = 300)

   # Distribution plot
   plot_distribution(cohort$SE, main = paste('Standard Errors of', prefixes[i]), xlab = 'SE')
   plot_distribution(cohort$P_VAL, main = paste('p-values of', prefixes[i]), xlab = 'p-value')

   # QQ plot.
   qqman::qq(cohort$P_VAL, main = sprintf('QQ plot of %s (lambda = %f)', prefixes[i], lambda=get_lambda(cohort,"P_VAL")))

   # Volcano plot.
   plot_volcano(cohort, "BETA", "P_VAL", main =paste('Volcano plot of', prefixes[i]) )

   dev.off()

   # Add mandatory data for precissionplot
   medianSE[i] <-  median(cohort$SE)
   value_N[i] <- N[i]
   cohort_label[i] <-  prefixes[i]

   # if n is defined for dichotomic condition :
   if(length(n) == length(N))  value_n[i] <- n[i]
}


# Data for Precission Plot
precplot.data <- as.data.frame(cbind( SE = medianSE, invSE = (1/medianSE), N = value_N, sqrt_N = sqrt(N), cohort = cohort_label ))

if(length(n) == length(N))
   precplot.data.n <- as.data.frame(cbind( SE = medianSE, invSE = (1/medianSE), N = value_n, sqrt_N = sqrt(n), cohort = cohort_label ))



##  Post model analysis  ##

if ( length(files) > 1)
{
   # Precission_Plot(N)
   plot_precissionp(precplot.data, paste(results_folder,"precision_SE_N.png",sep='/'), main = "Precision Plot - 1/median(SE) vs sqrt(N)")

   # Precission_Plot(n)
   if(length(n) == length(N))
      plot_precissionp(precplot.data.n, paste(results_folder,"precision_SE_N.png"), main = "Subgroup Precision Plot -  1/median(SE) vs sqrt(n)")

   # Venn_Diagrams()
   for (i in 1:length(venn_diagrams))
      plot_venndiagram(venn_diagrams[[i]], qcpath = results_folder, plotpath = results_folder )

}


## ############### ##
##  Meta-Analysis  ##
## ############### ##

## -- Variable definition -- ##

# Define data for each meta-analysis
metafiles <- list(
   'MetaA1' = c('PACE_AQUA_A1','PACE_IMMA_A1', 'RICHS_A1' ),
   'MetaA2' = c('PACE_AQUA_A2','PACE_IMMA_A2', 'RICHS_A2' ),
   'MetaB' = c('PACE_IMMA_B1','PACE_IMMA_B2'))

# Define maximum percent missing for each CpG
pcentMissing <- 0.8 # CpGs with precense lower than pcentMissing after GWAS meta-analysis will be deleted from the study.

folder <- results_folder

## Create directory for GWAMA configuration files and GWAMA_Results
if(!dir.exists(file.path(getwd(), paste(results_folder, "GWAMA", sep="/") )))
   suppressWarnings(dir.create(file.path(getwd(), paste(results_folder, "GWAMA", sep="/"))))

## Create directory for GWAMA_Results
outputfolder <- paste0(results_folder, "/GWAMA_Results")
if(!dir.exists(file.path(getwd(), outputfolder )))
   suppressWarnings(dir.create(file.path(getwd(), outputfolder)))


# Create map file for GWAMA --> Used in Manhattan plots
hapmapfile <- paste(results_folder,"GWAMA", "hapmap.map" ,sep = "/")
generate_hapmap_file(artype, hapmapfile)

# GWAMA binary path
#.Original.# gwama.dir <- paste0(Sys.getenv("HOME"), "/data/EWAS_metaanalysis/1_QC_results_cohorts/GWAMA/")
gwama.dir <- "/Users/mailos/tmp/GWAMA_v2/"

for( metf in 1:length(metafiles))
{

   list.lowCpGs <- NULL

   # Create folder for a meta-analysis in GWAMA folder, here we store the GWAMA input files for each meta-analysis,
   # We create one for complete meta-analysis
   if(!dir.exists(file.path(getwd(), paste(results_folder,"GWAMA", names(metafiles)[metf] ,sep="/") )))
      suppressWarnings(dir.create(file.path(getwd(), paste(results_folder,"GWAMA", names(metafiles)[metf], sep="/"))))
   # We create another for meta-analysis without filtered CpGs with low percentage (sufix _Filtr)
   if(!dir.exists(file.path(getwd(), paste0(results_folder,"/GWAMA/", names(metafiles)[metf],"_Filtr") )))
      suppressWarnings(dir.create(file.path(getwd(), paste0(results_folder,"/GWAMA/", names(metafiles)[metf],"_Filtr"))))

   # GWAMA File name base
   inputfolder <- paste0(results_folder,"/GWAMA/",  names(metafiles)[metf])

   modelfiles <- unlist(metafiles[metf])

   runs <- c('Normal', 'lowcpgs') # Execution with all CpGs and without filtered CpGs
   lowCpGs = FALSE;
   outputfiles <- list()

   outputgwama <- paste(outputfolder,names(metafiles)[metf],sep = '/')

   for(j in 1:length(runs))
   {
      if(runs[j]=='lowcpgs') {
         lowCpGs = TRUE
         # Get low presence CpGs in order to exclude this from the new meta-analysis
         list.lowCpGs <- get_low_presence_CpGs(outputfiles[[j-1]], pcentMissing)
         inputfolder <- paste0(results_folder,"/GWAMA/",  names(metafiles)[metf], "_Filtr")
         outputgwama <- paste0(outputgwama,"_Filtr")
      }

      # Create GWAMA files for each file in meta-analysis and execute GWAMA
      for ( i in 1:length(modelfiles) )
         create_GWAMA_files(results_folder,  modelfiles[i], inputfolder, N[i], list.lowCpGs )

      #.Original.#outputfiles[[runs[j]]] <- execute_GWAMA_MetaAnalysis(prefixgwama, names(metafiles)[metf])
      outputfiles[[runs[j]]] <- run_GWAMA_MetaAnalysis(inputfolder, outputgwama, names(metafiles)[metf], gwama.dir)

      # Post-metha-analysis QC --- >>> adds BN and FDR adjustment
      ##..## dataPost <- get_descriptives_postGWAMA(outputfolder, outputfiles[[runs[j]]], modelfiles, names(metafiles)[metf], artype )
      dataPost <- get_descriptives_postGWAMA(outputgwama, outputfiles[[runs[j]]], modelfiles, names(metafiles)[metf], artype )

      # Forest-Plot
      plot_ForestPlot( dataPost, metafiles[[metf]], runs[j], inputfolder, names(metafiles)[metf]  )

   }


}




## ############ ##
##  Enrichment  ##
## ############ ##



FilesToEnrich <- c('./QC_Results/toenrich/CpGstoEnrich.txt', # new line separation
                   './QC_Results/GWAMA_Results/MetaA1/MetaA1_Fixed_Modif.out' # file with GWAMA adjusted results
                   )
BN <-  TRUE
FDR <- 0.7
pvalue <- 0.05


## Check if we have any files to enrich and if these files exists
if (length(FilesToEnrich)>=1 & FilesToEnrich[1]!='') {
   for ( i in 1:length(FilesToEnrich))
      if (!file.exists(FilesToEnrich[i])) stop(paste0('File ',FilesToEnrich[i],' does not exsits, please check file' ))
}

outputfolder <- file.path(getwd(), "Enrichment" )


# Create dir to put results from enrichment
if(!dir.exists(outputfolder ))
   suppressWarnings(dir.create(outputfolder))

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
      if(dim(data)[1] <= 1 | dim(data)[2] <= 1) {
         data <- read.table(FilesToEnrich[i], dec = ".") # Avoid header
         data <- as.vector(t(data))
         data <- get_annotattions(data, artype, FilesToEnrich[i], outputfolder )
         allCpGs <- TRUE
      }

      ## -- Functional Enrichmnet
      ## ------------------------

      # Enrichment with missMethyl - GO and KEGG --> Writes results to outputfolder
      miss_enrich <- missMethyl_enrichment(data, outputfolder, FilesToEnrich[i], artype, BN, FDR, pvalue, allCpGs, plots = TRUE)



      ## -- Online Tools

      #     - Consensus path http://cpdb.molgen.mpg.de/ (gene-set analysis – over-representation analysis)
      #     - enrichr https://amp.pharm.mssm.edu/Enrichr/ (it alows gene-set analysis, but also disease
      #       enrichment and enrichemnt for transcription factors)



      ## -- Molecular Enrichmnet
      ## -----------------------

      # Molecular Signatures Database enrichment
      msd_enrich <- MSigDB_enrichment(data, outputfolder, FilesToEnrich[i], artype, BN, FDR, pvalue, allCpGs)


      # get unique genes from data
      geneUniv <- lapply( lapply(miss_enrich[grepl("signif", names(miss_enrich))], function(cpgs) { data[which(as.character(data$CpGs) %in% cpgs),]$UCSC_RefGene_Name}), getUniqueGenes)

      ## --  FER - with Online tools ==> A veure si es pot fer l'enriquiment automàtic amb les eines web.....
      ##### TO DO  : Fer-ho amb script automàtic??? --> Descarregar la informació associada als gens des de les webs utilitzant
      #####          els scripts que proporciona la web?? , mirar com s'han d'utilitzar i amb quin llenguatge.



      ## --  CpG Gene position - Fisher Test
      ## -----------------------------------

      get_descriptives_GenePosition(crom_data$UCSC_RefGene_Group, crom_data$Bonferroni, "Bonferroni", outputdir = "Fisher_BN_GenePosition_Desc", outputfile = FilesToEnrich[i])
      get_descriptives_GenePosition(crom_data$UCSC_RefGene_Group, crom_data$bFDR , "FDR", outputdir = "Fisher_FDR_GenePosition_Desc", outputfile = FilesToEnrich[i])

      # OR fir FDR significatives OR by position Relative to Island
      relative_island_fdr <- getAllFisherTest(crom_data$bFDR, crom_data$UCSC_RefGene_Group,
                                              outputdir = "Fisher_FDR_GenePosition", outputfile = FilesToEnrich[i], plots = TRUE )
      # OR fir FDR significatives OR by position Relative to Island
      relative_island_fdr_hyper <- getAllFisherTest(FDR_Hyper, crom_data$UCSC_RefGene_Group,
                                                    outputdir = "Fisher_FDRHyper_GenePosition", outputfile = FilesToEnrich[i], plots = TRUE )
      # OR fir FDR significatives OR by position Relative to Island
      relative_island_fdr_hypo <- getAllFisherTest(FDR_Hypo, crom_data$UCSC_RefGene_Group,
                                                   outputdir = "Fisher_FDRHypo_GenePosition", outputfile = FilesToEnrich[i], plots = TRUE )


      ## --  CpG Gene position - HyperGeometric Test
      ## -------------------------------------------

      # http://mengnote.blogspot.com.es/2012/12/calculate-correct-hypergeometric-p.html

      # Get hyper-geometric test for each Island relative position :
      # Depletion and Enrichment , Depletion and Enrichment for Hyper and Depletion and Enrichment for Hypo

      # Hypergeometric Test for FDR significatives by position Relative to Island
      hypergeo_relisland_fdr <- getAllHypergeometricTest(crom_data$bFDR, crom_data$UCSC_RefGene_Group, outputdir = "HyperG_FDR_GenePosition", outputfile = FilesToEnrich[i])
      # Hypergeometric Test for FDR Hyper significatives by position Relative to Island
      hypergeo_relisland_fdrhyper <- getAllHypergeometricTest(FDR_Hyper, crom_data$UCSC_RefGene_Group, outputdir = "HyperG_FDRHyper_GenePosition", outputfile = FilesToEnrich[i])
      # Hypergeometric Test for FDR Hypo significatives by position Relative to Island
      hypergeo_relisland_fdrhypo <- getAllHypergeometricTest(FDR_Hypo, crom_data$UCSC_RefGene_Group, outputdir = "HyperG_FDRHypo_GenePosition", outputfile = FilesToEnrich[i])





      ## --  CpG Island relative position - Fisher Test
      ## ----------------------------------------------

      get_descriptives_RelativetoIsland(crom_data$Relation_to_Island, crom_data$Bonferroni, "Bonferroni", outputdir = "Fisher_BN_RelativeToIsland", outputfile = FilesToEnrich[i])
      get_descriptives_RelativetoIsland(crom_data$Relation_to_Island, crom_data$bFDR , "FDR", outputdir = "Fisher_FDR_RelativeToIsland", outputfile = FilesToEnrich[i])

      # OR fir FDR significatives OR by position Relative to Island
      relative_island_fdr <- getAllFisherTest(crom_data$bFDR, crom_data$Relation_to_Island,
                                    outputdir = "Fisher_FDR_RelativeToIsland", outputfile = FilesToEnrich[i], plots = TRUE )
      # OR fir FDR significatives OR by position Relative to Island
      relative_island_fdr_hyper <- getAllFisherTest(FDR_Hyper, crom_data$Relation_to_Island,
                                                    outputdir = "Fisher_FDRHyper_RelativeToIsland", outputfile = FilesToEnrich[i], plots = TRUE )
      # OR fir FDR significatives OR by position Relative to Island
      relative_island_fdr_hypo <- getAllFisherTest(FDR_Hypo, crom_data$Relation_to_Island,
                                                    outputdir = "Fisher_FDRHypo_RelativeToIsland", outputfile = FilesToEnrich[i], plots = TRUE )


      ## --  CpG Island relative position - HyperGeometric Test
      ## ------------------------------------------------------
      # http://mengnote.blogspot.com.es/2012/12/calculate-correct-hypergeometric-p.html

      # Get hyper-geometric test for each Island relative position :
      # Depletion and Enrichment , Depletion and Enrichment for Hyper and Depletion and Enrichment for Hypo

      # Hypergeometric Test for FDR significatives by position Relative to Island
      hypergeo_relisland_fdr <- getAllHypergeometricTest(crom_data$bFDR, crom_data$Relation_to_Island, outputdir = "HyperG_FDR_RelativeToIsland", outputfile = FilesToEnrich[i])
      # Hypergeometric Test for FDR Hyper significatives by position Relative to Island
      hypergeo_relisland_fdrhyper <- getAllHypergeometricTest(FDR_Hyper, crom_data$Relation_to_Island, outputdir = "HyperG_FDRHyper_RelativeToIsland", outputfile = FilesToEnrich[i])
      # Hypergeometric Test for FDR Hypo significatives by position Relative to Island
      hypergeo_relisland_fdrhypo <- getAllHypergeometricTest(FDR_Hypo, crom_data$Relation_to_Island, outputdir = "HyperG_FDRHypo_RelativeToIsland", outputfile = FilesToEnrich[i])



      ## --  ROADMAP  -  Metilation in Cromatine States - BLOOD
      ## -------------------------------------------------------

      ##       Analysis of methylation changes in the different chromatin states (CpGs are diff meth in some states and others don't)

      # Prepare data
      crom_data <- addCrom15Columns(data, "rs_number") # Adds chromatine state columns
      crom_data$meth_state <- getHyperHypo(data$beta) # Classify methylation into Hyper and Hypo
      crom_data$bFDR <- getBinaryClassificationYesNo(crom_data$FDR, "<", FDR) # Classify fdr into "yes" and no taking into account FDR significance level

      # Columns with chromatin status information :
      ChrStatCols <- c("TssA","TssAFlnk","TxFlnk","TxWk","Tx","EnhG","Enh","ZNF.Rpts","Het","TssBiv","BivFlnk","EnhBiv","ReprPC","ReprPCWk","Quies")

      ## -- Define FDR and BN filter with Hyper and Hypo data

      # CpGs FDR and Hyper and Hypo respectively
      FDR_Hyper <- ifelse(crom_data$bFDR == 'yes' & crom_data$meth_state=='Hyper', "yes", "no")
      FDR_Hypo <- ifelse(crom_data$bFDR == 'yes' & crom_data$meth_state=='Hypo', "yes", "no")

      # CpGs Bonferroni and Hyper and Hypo respectively
      BN_Hyper <- ifelse(crom_data$Bonferroni == 'yes' & crom_data$meth_state=='Hyper', "yes", "no")
      BN_Hypo <- ifelse(crom_data$Bonferroni == 'yes' & crom_data$meth_state=='Hypo', "yes", "no")



      # FDR significatives regression by chromatin state
      chrom_states_fdr <- getAllChromStateOR(crom_data$bFDR, crom_data[,ChrStatCols],outputdir = "OR_FDR_States", outputfile = FilesToEnrich[i], plots = TRUE )
      # FDR significative with Hypermetilation vs no significative hypomethylation
      chrom_states_fdr_hyper <- getAllChromStateOR(FDR_Hyper, crom_data[,ChrStatCols], outputdir = "OR_FDRHyper_States", outputfile = FilesToEnrich[i], plots = TRUE )
      # FDR significative with Hypo-methylation vs no significative hyper-methylation
      chrom_states_fdr_hypo <- getAllChromStateOR(FDR_Hypo, crom_data[,ChrStatCols], outputdir = "OR_FDRHypo_States", outputfile = FilesToEnrich[i], plots = TRUE )



      ## -- ROADMAP  -  Regulatory feature enrichment analysis - PLACENTA
      ## -----------------------------------------------------------------

      # Convert to Genomic Ranges
      data.GRange <- GRanges(
         seqnames = Rle(data$chr),
         ranges=IRanges(data$pos, end=data$pos),
         name=data$rs_number,
         chr=data$chromosome,
         pos=data$pos
         )
      names(data.GRange) <- data.GRange$name

      # Find overlaps between CpGs and Fetal Placenta (States 15 ans 18)
      over15 <- findOverlapValues(data.GRange, FP_15_E091 )
      over18 <- findOverlapValues(data.GRange, FP_18_E091 )

      # Add states 15 and 18 to data.GRange file and write to a file : CpGs, state15 and state18
      data.chrstates <- c(mcols(over15$ranges), over15$values, over18$values)
      colnames(data.chrstates)[grep("States",colnames(data.chrstates))] <-  c("States15_FP", "States18_FP")
      # Merge annotated data with chromatine states with states with data
      crom_data <- merge(crom_data, data.chrstates, by.x = "rs_number", by.y = "name" )

      fname <- paste0("Regulatory_feature_enrichment/List_CpGs_",
                      tools::file_path_sans_ext(basename(FilesToEnrich[i])),
                      "_annot_plac_chr_states.txt")
      dir.create("Regulatory_feature_enrichment", showWarnings = FALSE)
      write.table( crom_data, fname, quote=F, row.names=F, sep="\t")

      ## --  HyperGeometric Test - States15_FP - BN

      # Hypergeometric Test for BN significatives by States15_FP
      hypergeo_States15FP_bn <- getAllHypergeometricTest(crom_data$Bonferroni, crom_data$States15_FP, outputdir = "HyperG_BN_States15_FP", outputfile = FilesToEnrich[i])
      # Hypergeometric Test for FDR significatives by States15_FP
      hypergeo_States15FP_bnhyper <- getAllHypergeometricTest(BN_Hyper, crom_data$States15_FP, outputdir = "HyperG_BNHyper_States15_FP", outputfile = FilesToEnrich[i])
      # Hypergeometric Test for FDR significatives by States15_FP
      hypergeo_States15FP_bnhypo <- getAllHypergeometricTest(BN_Hypo, crom_data$States15_FP, outputdir = "HyperG_BNHypo_States15_FP", outputfile = FilesToEnrich[i])

      ## --  Resume in a table - HyperGeometric Test - States15_FP - BN
      resdata <- summary_States_FP_Table( crom_data$Bonferroni, BN_Hyper, BN_Hypo, crom_data$States15_FP, outputdir = "Sum_States15_FP", outputfile = FilesToEnrich[i], plot = TRUE )


      ## --  HyperGeometric Test - States18_FP - BN

      # Hypergeometric Test for BN significatives by States15_FP
      hypergeo_States15FP_bn <- getAllHypergeometricTest(crom_data$Bonferroni, crom_data$States18_FP, outputdir = "HyperG_BN_States18_FP", outputfile = FilesToEnrich[i])
      # Hypergeometric Test for FDR significatives by States15_FP
      hypergeo_States15FP_bnhyper <- getAllHypergeometricTest(BN_Hyper, crom_data$States18_FP, outputdir = "HyperG_BNHyper_States18_FP", outputfile = FilesToEnrich[i])
      # Hypergeometric Test for FDR significatives by States15_FP
      hypergeo_States15FP_bnhypo <- getAllHypergeometricTest(BN_Hypo, crom_data$States18_FP, outputdir = "HyperG_BNHypo_States18_FP", outputfile = FilesToEnrich[i])

      ## --  Resume in a table - HyperGeometric Test - States15_FP - BN
      resdata <- summary_States_FP_Table( crom_data$Bonferroni, BN_Hyper, BN_Hypo, crom_data$States18_FP, outputdir = "Sum_States18_FP", outputfile = FilesToEnrich[i], plot = TRUE )



      ## -- Partially Methylated Domains (PMDs) PLACENTA
      ## ------------------------------------------------

      # Create genomic ranges from PMD data
      PMD.GRange <- getPMDGenomicRanges(PMD_placenta$Chr_PMD, PMD_placenta$Start_PMD, PMD_placenta$End_PMD)

      # Find overlaps between CpGs and PMD (find subject hits, query hits )
      overPMD <- findOverlapValues(data.GRange, PMD.GRange )

      #Create a data.frame with CpGs and PMDs information
      mdata <- as.data.frame(cbind(DataFrame(CpG = data.GRange$name[overPMD$qhits]), DataFrame(PMD = PMD.GRange$name[overPMD$shits])))

      # Merge with results from meta-analysis (A2)
      crom_data <- merge(crom_data, mdata, by.x="rs_number", by.y="CpG",all=T)
      crom_data <- crom_data[order(crom_data$p.value),]


      ## --  HyperGeometric Test - PMD - BN

      # Hypergeometric Test for BN significatives by PMD
      hypergeo_States15FP_bn <- getAllHypergeometricTest(crom_data$Bonferroni, crom_data$PMD, outputdir = "HyperG_BN_PMD", outputfile = FilesToEnrich[i])
      # Hypergeometric Test for FDR significatives by PMD
      hypergeo_States15FP_bnhyper <- getAllHypergeometricTest(BN_Hyper, crom_data$PMD, outputdir = "HyperG_BNHyper_PMD", outputfile = FilesToEnrich[i])
      # Hypergeometric Test for FDR significatives by PMD
      hypergeo_States15FP_bnhypo <- getAllHypergeometricTest(BN_Hypo, crom_data$PMD, outputdir = "HyperG_BNHypo_PMD", outputfile = FilesToEnrich[i])

      ## --  HyperGeometric Test - PMD - FDR

      # Hypergeometric Test for FDR significatives by PMD
      hypergeo_States15FP_bn <- getAllHypergeometricTest(crom_data$bFDR, crom_data$PMD, outputdir = "HyperG_FDR_PMD", outputfile = FilesToEnrich[i])
      # Hypergeometric Test for FDR significatives by PMD
      hypergeo_States15FP_bnhyper <- getAllHypergeometricTest(FDR_Hyper, crom_data$PMD, outputdir = "HyperG_FDRHyper_PMD", outputfile = FilesToEnrich[i])
      # Hypergeometric Test for FDR significatives by PMD
      hypergeo_States15FP_bnhypo <- getAllHypergeometricTest(FDR_Hypo, crom_data$PMD, outputdir = "HyperG_FDRHypo_PMD", outputfile = FilesToEnrich[i])


      ## --  Filter CpGs in Islands, Shores and Promoters





      #
   }

} else{
   print ("Error no data to enrich.")
}
