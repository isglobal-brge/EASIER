## ############### ##
##  Meta-Analysis  ##
## ############### ##


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

# Install methyTools (if needed)
# devtools::install_github("isglobal-brge/EASIER@HEAD")
library(methyTools)


########## ----------  VARIABLES DEFINED BY USER  ----------  ##########

## -- Variable definition for Meta-Analysis -- ##

# Files used in QC, needed in meta-analysis to plot ForestPlot
files <- c('data/Cohort1_Model1_20170713.txt',
           'data/Cohort1_Model2_20170713.txt',
           'data/PROJ1_Cohort3_Model1_date_v2.txt',
           'data/PROJ1_Cohort3_Model2_date_v2.txt',
           'data/PROJ1_Cohort2_Plate_ModelA1_20170309.txt',
           'data/PROJ1_Cohort2_Plate_ModelA2_20170309.txt',
           'data/PROJ1_Cohort2_Plate_ModelB1_20170320.txt',
           'data/PROJ1_Cohort2_Plate_ModelB2_20170320.txt',
           'data/PROJ1_Cohort2_Plate_ModelC1_20170818.txt',
           'data/PROJ1_Cohort2_Plate_ModelC2_20170818.txt')

# Prefixes for each file
prefixes <- c('Cohort1_A1', 'Cohort1_A2',
              'PROJ1_Cohort2_A1','PROJ1_Cohort2_A2', 'PROJ1_Cohort2_B1', 'PROJ1_Cohort2_B2', 'PROJ1_Cohort2_C1', 'PROJ1_Cohort2_C2',
              'PROJ1_Cohort3_A1', 'P1_Cohort3_A2')

# Samples in original files used in QC
N <- c(100, 100, 166, 166, 166, 166, 166, 166, 240, 240 )

# Array type, used : EPIC or 450K
artype <- '450K'

# Define data for each meta-analysis
metafiles <- list(
   'MetaA1' = c('Cohort1_A1','PROJ1_Cohort2_A1', 'PROJ1_Cohort3_A1' ),
   'MetaA2' = c('Cohort1_A2','PROJ1_Cohort2_A2', 'P1_Cohort3_A2' ),
   'MetaB' = c('PROJ1_Cohort2_B1','PROJ1_Cohort2_B2'))

# Define maximum percent missing for each CpG
pcentMissing <- 0.8 # CpGs with precense lower than pcentMissing after GWAS meta-analysis will be deleted from the study.


# Paths with QCResults and path to store GWAMA results
results_folder <- 'QC_Results'
results_gwama <- '.'


# Venn diagrams
venn_diagrams <- list(
   c("MetaA1", "MetaA2", "MetaB" ),
   c("MetaA1_Filtr", "MetaA2_Filtr", "MetaB_Filtr" )
)

########## ----------  END VARIABLES DEFINED BY USER  ----------  ##########



# GWAMA binary path  (GWAMA IsGlobal Server installation)
gwama.dir <- paste0(Sys.getenv("HOME"), "/data/EWAS_metaanalysis/1_QC_results_cohorts/GWAMA/")
gwama.dir <- "/Users/mailos/tmp/GWAMA_v2/"

## Create directory for GWAMA configuration files and GWAMA_Results
if(!dir.exists(file.path(paste(results_gwama, "GWAMA", sep="/") )))
   suppressWarnings(dir.create(file.path( paste(results_gwama, "GWAMA", sep="/"))))

## Create directory for GWAMA_Results
outputfolder <- paste0(results_gwama, "/GWAMA_Results")
if(!dir.exists(file.path( outputfolder )))
   suppressWarnings(dir.create(file.path( outputfolder)))



# Create map file for GWAMA --> Used in Manhattan plots
hapmapfile <- paste(results_gwama,"GWAMA", "hapmap.map" ,sep = "/")
generate_hapmap_file(artype, hapmapfile)



for( metf in 1:length(metafiles))
{

   list.lowCpGs <- NULL

   # Create folder for a meta-analysis in GWAMA folder, here we store the GWAMA input files for each meta-analysis,
   # We create one for complete meta-analysis
   if(!dir.exists(file.path( paste(results_gwama,"GWAMA", names(metafiles)[metf] ,sep="/") )))
      suppressWarnings(dir.create(file.path( paste(results_gwama,"GWAMA", names(metafiles)[metf], sep="/"))))
   # We create another for meta-analysis without filtered CpGs with low percentage (sufix _Filtr)
   if(!dir.exists(file.path( paste0(results_gwama,"/GWAMA/", names(metafiles)[metf],"_Filtr") )))
      suppressWarnings(dir.create(file.path( paste0(results_gwama,"/GWAMA/", names(metafiles)[metf],"_Filtr"))))

   # GWAMA File name base
   inputfolder <- paste0(results_gwama,"/GWAMA/",  names(metafiles)[metf])

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
         inputfolder <- paste0(results_gwama,"/GWAMA/",  names(metafiles)[metf], "_Filtr")
         outputgwama <- paste0(outputgwama,"_Filtr")
      }

      # Create GWAMA files for each file in meta-analysis and execute GWAMA
      for ( i in 1:length(modelfiles) )
         create_GWAMA_files(file.path(results_folder,modelfiles[i]),  modelfiles[i], inputfolder, N[i], list.lowCpGs )

      #.Original.#outputfiles[[runs[j]]] <- execute_GWAMA_MetaAnalysis(prefixgwama, names(metafiles)[metf])
      outputfiles[[runs[j]]] <- run_GWAMA_MetaAnalysis(inputfolder, outputgwama, names(metafiles)[metf], gwama.dir)

      # Post-metha-analysis QC --- >>> adds BN and FDR adjustment
      ##..## dataPost <- get_descriptives_postGWAMA(outputfolder, outputfiles[[runs[j]]], modelfiles, names(metafiles)[metf], artype )
      dataPost <- get_descriptives_postGWAMA(outputgwama, outputfiles[[runs[j]]], modelfiles, names(metafiles)[metf], artype, N[which(prefixes %in% modelfiles)] )

      # Forest-Plot
      plot_ForestPlot( dataPost, metafiles[[metf]], runs[j], inputfolder, names(metafiles)[metf], files, outputgwama  )

   }

}


# Venn_Diagrams for for meta-analysis with fixed effects
for (i in 1:length(venn_diagrams))
   plot_venndiagram(venn_diagrams[[i]], qcpath = outputfolder, plotpath =  paste0(results_gwama, "/GWAMA_Results"), pattern = '_Fixed_Modif.out',bn='Bonferroni', fdr='FDR')


if(dir.exists(file.path( paste(results_gwama, "GWAMA", sep="/") )))
   unlink(file.path(results_gwama, "GWAMA"), recursive=TRUE)


