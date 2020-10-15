if (!require(rasterpdf, quietly = TRUE)) install.packages('rasterpdf')
if (!require(ggplot2, quietly = TRUE)) install.packages('ggplot2')
if (!require(VennDiagram, quietly = TRUE)) install.packages('VennDiagram')
if (!require(RColorBrewer, quietly = TRUE)) install.packages('RColorBrewer')
if (!require(tibble, quietly = TRUE)) install.packages('tibble')
if (!require(dplyr, quietly = TRUE)) install.packages('dplyr')
if (!require(stringr, quietly = TRUE)) install.packages('stringr')
if (!require(meta, quietly = TRUE)) install.packages('meta') # Forest Plot
if (!require(missMethyl, quietly = TRUE)) BiocManager::install( "missMethyl" )
if (!require(org.Hs.eg.db, quietly = TRUE)) BiocManager::install( "org.Hs.eg.db" )

if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")

#BiocManager::install( c("IlluminaHumanMethylation450kanno.ilmn12.hg19",
#                        "IlluminaHumanMethylation450kanno.ilmn12.hg19",
#                        "missMethyl",
#                        "org.Hs.eg.db") )

library(org.Hs.eg.db)

library(methyTools)

setwd("")


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


## Create directory for GWAMA configuration files and GWAMA_Results
if(!dir.exists(file.path(getwd(), paste(folder, "GWAMA", sep="/") )))
   suppressWarnings(dir.create(file.path(getwd(), paste(folder, "GWAMA", sep="/"))))

## Create directory for GWAMA_Results
outputfolder <- paste0(folder, "/GWAMA_Results")
if(!dir.exists(file.path(getwd(), outputfolder )))
   suppressWarnings(dir.create(file.path(getwd(), outputfolder)))


# Create map file for GWAMA --> Used in Manhattan plots
hapmapfile <- paste(folder,"GWAMA", "hapmap.map" ,sep = "/")
generate_hapmap_file(artype, hapmapfile)

# GWAMA binary path
#.Original.# gwama.dir <- paste0(Sys.getenv("HOME"), "/data/EWAS_metaanalysis/1_QC_results_cohorts/GWAMA/")
gwama.dir <- "/Users/mailos/tmp/GWAMA_v2/"

for( metf in 1:length(metafiles))
{

   list.lowCpGs <- NULL

   # Create folder for a meta-analysis in GWAMA folder, here we store the GWAMA input files for each meta-analysis,
   # We create one for complete meta-analysis
   if(!dir.exists(file.path(getwd(), paste(folder,"GWAMA", names(metafiles)[metf] ,sep="/") )))
      suppressWarnings(dir.create(file.path(getwd(), paste(folder,"GWAMA", names(metafiles)[metf], sep="/"))))
   # We create another for meta-analysis without filtered CpGs with low percentage (sufix _Filtr)
   if(!dir.exists(file.path(getwd(), paste0(folder,"/GWAMA/", names(metafiles)[metf],"_Filtr") )))
      suppressWarnings(dir.create(file.path(getwd(), paste0(folder,"/GWAMA/", names(metafiles)[metf],"_Filtr"))))

   # GWAMA File name base
   inputfolder <- paste0(folder,"/GWAMA/",  names(metafiles)[metf])

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
         inputfolder <- paste0(folder,"/GWAMA/",  names(metafiles)[metf], "_Filtr")
         outputgwama <- paste0(outputgwama,"_Filtr")
      }

      # Create GWAMA files for each file in meta-analysis and execute GWAMA
      for ( i in 1:length(modelfiles) )
         create_GWAMA_files(folder,  modelfiles[i], inputfolder, N[i], list.lowCpGs )

      #.Original.#outputfiles[[runs[j]]] <- execute_GWAMA_MetaAnalysis(prefixgwama, names(metafiles)[metf])
      outputfiles[[runs[j]]] <- run_GWAMA_MetaAnalysis(inputfolder, outputgwama, names(metafiles)[metf], gwama.dir)

      # Post-metha-analysis QC --- >>> adds BN and FDR adjustment
      ##..## dataPost <- get_descriptives_postGWAMA(outputfolder, outputfiles[[runs[j]]], modelfiles, names(metafiles)[metf], artype )
      dataPost <- get_descriptives_postGWAMA(outputgwama, outputfiles[[runs[j]]], modelfiles, names(metafiles)[metf], artype )

      # Forest-Plot
      plot_ForestPlot( dataPost, metafiles[[metf]], runs[j], inputfolder, names(metafiles)[metf]  )

   }


}

