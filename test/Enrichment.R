## ############ ##
##  Enrichment  ##
## ############ ##


# Install requirede libraries
if (!require(rasterpdf, quietly = TRUE)) install.packages('rasterpdf', repos = 'https://cran.rediris.es/' )
if (!require(meta, quietly = TRUE)) install.packages('meta', repos = 'https://cran.rediris.es/' )
if (!require(tibble, quietly = TRUE)) install.packages('tibble')
if (!require(dplyr, quietly = TRUE)) install.packages('dplyr')
if (!require(tidyverse, quietly = TRUE)) install.packages::install( "tidyverse" )
if (!require(stringr, quietly = TRUE)) install.packages('stringr')
if (!require(meta, quietly = TRUE)) install.packages('meta') # Forest Plot
if (!require(ggplot2, quietly = TRUE)) install.packages('ggplot2')
if (!require(VennDiagram, quietly = TRUE)) install.packages('VennDiagram')
if (!require(RColorBrewer, quietly = TRUE)) install.packages('RColorBrewer')
if (!require(reshape, quietly = TRUE)) install.packages('reshape')
if (!require(ggsignif, quietly = TRUE)) install.packages('ggsignif')
if (!require(tools, quietly = TRUE)) install.packages('tools')

# Bioconductor
if (!require(missMethyl, quietly = TRUE)) BiocManager::install( "missMethyl" )
if (!require(org.Hs.eg.db, quietly = TRUE)) BiocManager::install( "org.Hs.eg.db" )
if (!require(GenomicRanges, quietly = TRUE)) BiocManager::install( "GenomicRanges" )
if (!require(rtracklayer, quietly = TRUE)) BiocManager::install( "rtracklayer" )

# Load libraries

library(rasterpdf)
library(meta)
library(tibble)
library(dplyr)
library(stringr)
library(meta)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
library(reshape)
library(ggsignif)

# Bioconductor
library(missMethyl)
library(org.Hs.eg.db)
library(GenomicRanges)
library(rtracklayer)



# Install methyTools (if needed)
# devtools::install_github("isglobal-brge/methyTools@HEAD")
library(methyTools)


# setwd("/Users/mailos/tmp/proves")

FilesToEnrich <- c('./GWAMA_Results/MetaA1/MetaA1_Fixed_Modif.out',
                   'GWAMA_Results/MetaA2/MetaA2_Fixed_Modif.out',
                   'toenrich/CpGstoEnrich.txt'
                   )

# Values for adjustment
BN <-  TRUE       # Use Bonferroni ?
FDR <- 0.7        # significance level for adjustment, if NA FDR is not used
pvalue <- 0.05    # significance level for p-value, if NA p-value is not used


# Array type, used : EPIC or 450K
artype <- '450K'


## Check if we have any files to enrich and if these files exists
if (length(FilesToEnrich)>=1 & FilesToEnrich[1]!='') {
   for ( i in 1:length(FilesToEnrich))
      if (!file.exists(FilesToEnrich[i])) stop(paste0('File ',FilesToEnrich[i],' does not exsits, please check file' ))
}


# Result paths definition for QC, Meta-Analysis and Enrichment
results_folder <- 'QC_Results'
results_gwama <- '.'
results_enrich <- 'Enrichment'


# Convert relative paths to absolute paths for FilesToEnrich
FilesToEnrich <- unlist(sapply(FilesToEnrich, function(file) { if(substr(file,1,1)!='.' & substr(file,1,1)!='/') file <- paste0('./',file) else file }))
FilesToEnrich <- sapply(FilesToEnrich, file_path_as_absolute)

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
      if(dim(data)[1] <= 1 | dim(data)[2] <= 1) {
         data <- read.table(FilesToEnrich[i], dec = ".") # Avoid header
         data <- as.vector(t(data))
         data <- get_annotattions(data, artype, FilesToEnrich[i], outputfolder )
         allCpGs <- TRUE
         data$chromosome <- substr(data$chr,4,length(data$chr))
         data$rs_number <- data$CpGs
      }

      ## -- Functional Enrichmnet
      ## ------------------------

      # Enrichment with missMethyl - GO and KEGG --> Writes results to outputfolder
      miss_enrich <- missMethyl_enrichment(data, outputfolder, FilesToEnrich[i], artype, BN, FDR, pvalue, allCpGs, plots = TRUE )



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



      if("FDR" %in% colnames(data) & "Bonferroni" %in% colnames(data))
      {
         ## --  CpG Gene position - Fisher Test
         ## -----------------------------------

         # Add column bFDR to data for that CpGs that accomplish with FDR
         data$bFDR <- getBinaryClassificationYesNo(data$FDR, "<", FDR) # Classify fdr into "yes" and no taking into account FDR significance level

         # Classify by Hyper and Hypo methylated
         data$meth_state <- getHyperHypo(data$beta) # Classify methylation into Hyper and Hypo

         # CpGs FDR and Hyper and Hypo respectively
         FDR_Hyper <- ifelse(data$bFDR == 'yes' & data$meth_state=='Hyper', "yes", "no")
         FDR_Hypo <- ifelse(data$bFDR == 'yes' & data$meth_state=='Hypo', "yes", "no")

         # CpGs Bonferroni and Hyper and Hypo respectively
         BN_Hyper <- ifelse(data$Bonferroni == 'yes' & data$meth_state=='Hyper', "yes", "no")
         BN_Hypo <- ifelse(data$Bonferroni == 'yes' & data$meth_state=='Hypo', "yes", "no")

         # Get descriptives
         get_descriptives_GenePosition(data$UCSC_RefGene_Group, data$Bonferroni, "Bonferroni", outputdir = "GenePosition/Fisher_BN_Desc", outputfile = FilesToEnrich[i])
         get_descriptives_GenePosition(data$UCSC_RefGene_Group, data$bFDR , "FDR", outputdir = "GenePosition/Fisher_FDR_Desc", outputfile = FilesToEnrich[i])

         ## --  Fisher Test - Gene position - FDR, FDR_hyper and FDR_hypo
         GenePosition_fdr <- getAllFisherTest(data$bFDR, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_FDR", outputfile = FilesToEnrich[i], plots = TRUE )
         GenePosition_fdr_hyper <- getAllFisherTest(FDR_Hyper, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_FDRHyper", outputfile = FilesToEnrich[i], plots = TRUE )
         GenePosition_fdr_hypo <- getAllFisherTest(FDR_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/Fisher_FDRHypo", outputfile = FilesToEnrich[i], plots = TRUE )


         ## --  CpG Gene position - HyperGeometric Test (Depletion and Enrichment)
         ## -------------------------------------------

         # http://mengnote.blogspot.com.es/2012/12/calculate-correct-hypergeometric-p.html

         # Get hyper-geometric test for each Island relative position :
         # Depletion and Enrichment , Depletion and Enrichment for Hyper and Depletion and Enrichment for Hypo

         ## --  HyperGeometric Test - Gene position - FDR, FDR_hyper and FDR_hypo

         hypergeo_relisland_fdr <- getAllHypergeometricTest(data$bFDR, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_FDR", outputfile = FilesToEnrich[i])
         hypergeo_relisland_fdrhyper <- getAllHypergeometricTest(FDR_Hyper, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_FDRHyper", outputfile = FilesToEnrich[i])
         hypergeo_relisland_fdrhypo <- getAllHypergeometricTest(FDR_Hypo, data$UCSC_RefGene_Group, outputdir = "GenePosition/HyperG_FDRHypo", outputfile = FilesToEnrich[i])





         ## --  CpG Island relative position - Fisher Test (Depletion and Enrichment)
         ## ----------------------------------------------

         get_descriptives_RelativetoIsland(data$Relation_to_Island, data$Bonferroni, "Bonferroni", outputdir = "RelativeToIsland/Fisher_BN_RelativeToIsland", outputfile = FilesToEnrich[i])
         get_descriptives_RelativetoIsland(data$Relation_to_Island, data$bFDR , "FDR", outputdir = "RelativeToIsland/Fisher_FDR_RelativeToIsland", outputfile = FilesToEnrich[i])


         ## --  Fisher Test - Position Relative to Island - FDR, FDR_hyper and FDR_hypo

         relative_island_fdr <- getAllFisherTest(data$bFDR, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_FDR", outputfile = FilesToEnrich[i], plots = TRUE )
         relative_island_fdr_hyper <- getAllFisherTest(FDR_Hyper, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_FDRHyper", outputfile = FilesToEnrich[i], plots = TRUE )
         relative_island_fdr_hypo <- getAllFisherTest(FDR_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/Fisher_FDRHypo", outputfile = FilesToEnrich[i], plots = TRUE )


         ## --  CpG Island relative position - HyperGeometric Test (Depletion and Enrichment)
         ## ------------------------------------------------------
         # http://mengnote.blogspot.com.es/2012/12/calculate-correct-hypergeometric-p.html

         # Get hyper-geometric test for each Island relative position :
         # Depletion and Enrichment , Depletion and Enrichment for Hyper and Depletion and Enrichment for Hypo

         ## --  HyperGeometric Test - Position Relative to Island - FDR, FDR_hyper and FDR_hypo

         hypergeo_relisland_fdr <- getAllHypergeometricTest(data$bFDR, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_FDR", outputfile = FilesToEnrich[i])
         hypergeo_relisland_fdrhyper <- getAllHypergeometricTest(FDR_Hyper, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_FDRHyper", outputfile = FilesToEnrich[i])
         hypergeo_relisland_fdrhypo <- getAllHypergeometricTest(FDR_Hypo, data$Relation_to_Island, outputdir = "RelativeToIsland/HyperG_FDRHypo", outputfile = FilesToEnrich[i])

      }



      ## --  ROADMAP  -  Metilation in Cromatine States - BLOOD
      ## -------------------------------------------------------

      ##       Analysis of methylation changes in the different chromatin states (CpGs are diff meth in some states and others don't)

      # Prepare data
      crom_data <- addCrom15Columns(data, "rs_number") # Adds chromatine state columns
      # crom_data$meth_state <- getHyperHypo(data$beta) # Classify methylation into Hyper and Hypo

      if("FDR" %in% colnames(data) & "Bonferroni" %in% colnames(data))
      {

         # Columns with chromatin status information :
         ChrStatCols <- c("TssA","TssAFlnk","TxFlnk","TxWk","Tx","EnhG","Enh","ZNF.Rpts","Het","TssBiv","BivFlnk","EnhBiv","ReprPC","ReprPCWk","Quies")

         ## -- Define FDR and BN filter with Hyper and Hypo data (Depletion and Enrichment)

         # CpGs FDR and Hyper and Hypo respectively
         FDR_Hyper <- ifelse(crom_data$bFDR == 'yes' & crom_data$meth_state=='Hyper', "yes", "no")
         FDR_Hypo <- ifelse(crom_data$bFDR == 'yes' & crom_data$meth_state=='Hypo', "yes", "no")

         # CpGs Bonferroni and Hyper and Hypo respectively
         BN_Hyper <- ifelse(crom_data$Bonferroni == 'yes' & crom_data$meth_state=='Hyper', "yes", "no")
         BN_Hypo <- ifelse(crom_data$Bonferroni == 'yes' & crom_data$meth_state=='Hypo', "yes", "no")

         ## --  HyperGeometric Test - Chromatin State - FDR, FDR_hyper and FDR_hypo (Depletion and Enrichment)

         chrom_states_fdr <- getAllChromStateOR(crom_data$bFDR, crom_data[,ChrStatCols],outputdir = "CromStates/OR_FDR", outputfile = FilesToEnrich[i], plots = TRUE )
         chrom_states_fdr_hyper <- getAllChromStateOR(FDR_Hyper, crom_data[,ChrStatCols], outputdir = "CromStates/OR_FDRHyper", outputfile = FilesToEnrich[i], plots = TRUE )
         chrom_states_fdr_hypo <- getAllChromStateOR(FDR_Hypo, crom_data[,ChrStatCols], outputdir = "CromStates/OR_FDRHypo", outputfile = FilesToEnrich[i], plots = TRUE )
      }


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

      if("FDR" %in% colnames(data) & "Bonferroni" %in% colnames(data))
      {

         ## --  HyperGeometric Test - States18_FP - BN,  BN_hyper and BN_hypo (Depletion and Enrichment)

         hypergeo_States15FP_bn <- getAllHypergeometricTest(crom_data$Bonferroni, crom_data$States15_FP, outputdir = "States15_FP/HyperG_BN", outputfile = FilesToEnrich[i])
         hypergeo_States15FP_bnhyper <- getAllHypergeometricTest(BN_Hyper, crom_data$States15_FP, outputdir = "States15_FP/HyperG_BNHyper", outputfile = FilesToEnrich[i])
         hypergeo_States15FP_bnhypo <- getAllHypergeometricTest(BN_Hypo, crom_data$States15_FP, outputdir = "States15_FP/HyperG_BNHypo", outputfile = FilesToEnrich[i])

         ## --  Resume in a table - HyperGeometric Test - States15_FP - BN
         resdata <- summary_States_FP_Table( crom_data$Bonferroni, BN_Hyper, BN_Hypo, crom_data$States15_FP, outputdir = "States15_FP/Summary_HyperG_BN", outputfile = FilesToEnrich[i], plot = TRUE )


         ## --  HyperGeometric Test - States18_FP - BN,  BN_hyper and BN_hypo (Depletion and Enrichment)

         hypergeo_States15FP_bn <- getAllHypergeometricTest(crom_data$Bonferroni, crom_data$States18_FP, outputdir = "States18_FP/HyperG_BNP", outputfile = FilesToEnrich[i])
         hypergeo_States15FP_bnhyper <- getAllHypergeometricTest(BN_Hyper, crom_data$States18_FP, outputdir = "States18_FP/HyperG_BNHyper", outputfile = FilesToEnrich[i])
         hypergeo_States15FP_bnhypo <- getAllHypergeometricTest(BN_Hypo, crom_data$States18_FP, outputdir = "States18_FP/HyperG_BNHypo", outputfile = FilesToEnrich[i])

         ## --  Resume in a table - HyperGeometric Test - States15_FP - BN
         resdata <- summary_States_FP_Table( crom_data$Bonferroni, BN_Hyper, BN_Hypo, crom_data$States18_FP, outputdir = "States18_FP/Summary_HyperG_BN", outputfile = FilesToEnrich[i], plot = TRUE )
      }


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
      # crom_data <- crom_data[order(crom_data$p.value),]

      # CpGs with PMD as NA
      PMD_NaN <- ifelse(is.na(crom_data$PMD),'IsNA','NotNA' )

      if("FDR" %in% colnames(data) & "Bonferroni" %in% colnames(data))
      {

         ## --  HyperGeometric Test - PMD - BN,  BN_hyper and BN_hypo  (Full data ) (Depletion and Enrichment)

         hypergeo_PMD_bn <- getAllHypergeometricTest(crom_data$Bonferroni, PMD_NaN, outputdir = "PMD/HyperG_BN", outputfile = FilesToEnrich[i])
         hypergeo_PMD_bnhyper <- getAllHypergeometricTest(BN_Hyper, PMD_NaN, outputdir = "PMD/HyperG_BNHyper", outputfile = FilesToEnrich[i])
         hypergeo_PMD_bnhypo <- getAllHypergeometricTest(BN_Hypo, PMD_NaN, outputdir = "PMD/HyperG_BNHypo", outputfile = FilesToEnrich[i])
         # Summary
         resdata <- summary_HyperGeometrics_Table( crom_data$Bonferroni, BN_Hyper, BN_Hypo, PMD_NaN, outputdir = "PMD/Summary_HyperG_BN", outputfile = FilesToEnrich[i], plot = TRUE )

         ## --  HyperGeometric Test - PMD - FDR,  FDR_hyper and FDR_hypo  (Full data ) (Depletion and Enrichment)

         hypergeo_PMD_fdr <- getAllHypergeometricTest(crom_data$bFDR, PMD_NaN, outputdir = "PMD/HyperG_FDR", outputfile = FilesToEnrich[i])
         hypergeo_PMD_fdrhyper <- getAllHypergeometricTest(FDR_Hyper, PMD_NaN, outputdir = "PMD/HyperG_FDRHyper", outputfile = FilesToEnrich[i])
         hypergeo_PMD_fdrhypo <- getAllHypergeometricTest(FDR_Hypo, PMD_NaN, outputdir = "PMD/HyperG_FDRHypo", outputfile = FilesToEnrich[i])
         # Summary
         resdata <- summary_HyperGeometrics_Table( crom_data$bFDR, FDR_Hyper, FDR_Hypo, PMD_NaN, outputdir = "PMD/Summary_HyperG_FDR", outputfile = FilesToEnrich[i], plot = TRUE )


         ## --  Filter CpGs in Islands, Shores and Promoters

         fcrom_data <- crom_data[crom_data$rs_number %in% filterCpGs(crom_data$rs_number, crom_data$Relation_to_Island, c("Island", "N_Shore", "N_Shore")),]
         fcrom_data <- fcrom_data[fcrom_data$rs_number %in% filterCpGs(fcrom_data$rs_number, fcrom_data$UCSC_RefGene_Group, c("TSS200", "TSS1500")),]

         ## -- Define FDR and BN filter with Hyper and Hypo data for filtered CpGs

         # CpGs FDR, Hyper and Hypo respectively
         FDR_Hyper_f <- ifelse(fcrom_data$bFDR == 'yes' & fcrom_data$meth_state=='Hyper', "yes", "no")
         FDR_Hypo_f <- ifelse(fcrom_data$bFDR == 'yes' & fcrom_data$meth_state=='Hypo', "yes", "no")

         # CpGs Bonferroni, Hyper and Hypo respectively
         BN_Hyper_f <- ifelse(fcrom_data$Bonferroni == 'yes' & fcrom_data$meth_state=='Hyper', "yes", "no")
         BN_Hypo_f <- ifelse(fcrom_data$Bonferroni == 'yes' & fcrom_data$meth_state=='Hypo', "yes", "no")

         # CpGs with PMD as NA
         PMD_NaN_f <- ifelse(is.na(fcrom_data$PMD),'IsNA','NotNA' )

         ## --  HyperGeometric Test - PMD - BN,  BN_hyper and BN_hypo  (Filtered data ) (Depletion and Enrichment)

         hypergeo_PMD_bn_filt <- getAllHypergeometricTest(fcrom_data$Bonferroni, PMD_NaN_f, outputdir = "PMD/HyperG_BN_filtered", outputfile = FilesToEnrich[i])
         hypergeo_PMD_bnhyper_filt <- getAllHypergeometricTest(BN_Hyper_f, PMD_NaN_f, outputdir = "PMD/HyperG_BNHyper_filtered", outputfile = FilesToEnrich[i])
         hypergeo_PMD_bnhypo_filt <- getAllHypergeometricTest(BN_Hypo_f, PMD_NaN_f, outputdir = "PMD/HyperG_BNHypo_filtered", outputfile = FilesToEnrich[i])
         # Summary
         resdata <- summary_HyperGeometrics_Table( fcrom_data$Bonferroni, BN_Hyper_f, BN_Hypo_f, PMD_NaN_f, outputdir = "PMD/Summary_HyperG_BN_filtered", outputfile = FilesToEnrich[i], plot = TRUE )

         ## --  HyperGeometric Test - PMD - FDR,  FDR_hyper and FDR_hypo  (Filtered data ) (Depletion and Enrichment)

         hypergeo_PMD_fdr_filt <- getAllHypergeometricTest(fcrom_data$bFDR, PMD_NaN_f, outputdir = "PMD/HyperG_FDR_filtered", outputfile = FilesToEnrich[i])
         hypergeo_PMD_fdrhyper_filt <- getAllHypergeometricTest(FDR_Hyper_f, PMD_NaN_f, outputdir = "PMD/HyperG_FDRHyper_filtered", outputfile = FilesToEnrich[i])
         hypergeo_PMD_fdrhypo_filt <- getAllHypergeometricTest(FDR_Hypo_f, PMD_NaN_f, outputdir = "PMD/HyperG_FDRHypo_filtered", outputfile = FilesToEnrich[i])
         # Summary
         resdata <- summary_HyperGeometrics_Table( fcrom_data$bFDR, FDR_Hyper_f, FDR_Hypo_f, PMD_NaN_f, outputdir = "PMD/Summary_HyperG_FDR_filtered", outputfile = FilesToEnrich[i], plot = TRUE )

      }

      # WRITE FINAL ENRICHMENT DATA
      write.table( crom_data, paste0( getwd(), "/",tools::file_path_sans_ext(basename(FilesToEnrich[i])),"_Enriched.csv" ) , quote=F, row.names=F, sep="\t")
   }

} else{
   print ("Error no data to enrich.")
}
