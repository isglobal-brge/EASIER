if (!require(rasterpdf, quietly = TRUE)) install.packages('rasterpdf')
if (!require(ggplot2, quietly = TRUE)) install.packages('ggplot2')
if (!require(VennDiagram, quietly = TRUE)) install.packages('VennDiagram')
if (!require(RColorBrewer, quietly = TRUE)) install.packages('RColorBrewer')
if (!require(tibble, quietly = TRUE)) install.packages('tibble')
if (!require(dplyr, quietly = TRUE)) install.packages('dplyr')
if (!require(stringr, quietly = TRUE)) install.packages('stringr')
if (!require(meta, quietly = TRUE)) install.packages('meta') # Forest Plot
if (!require(reshape2, quietly = TRUE)) install.packages('reshape2') # Forest Plot

# devtools::install_github("isglobal-brge/EASIER@HEAD")

library(EASIER)

#..# setwd("~/Library/Mobile Documents/com~apple~CloudDocs/PROJECTES/Treballant/EASIER")
#..# setwd("/Users/mailos/tmp/proves")

########## ----------  VARIABLES DEFINED BY USER  ----------  ##########

# (this data is only an example)

# Set working directory to metaanalysis folder
setwd("<path to metaanalysis folder>/metaanalysis")

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

# Result folder
results_folder <- 'QC_Results'

# Prefixes for each file
prefixes <- c('Cohort1_A1', 'Cohort1_A2',
              'PROJ1_Cohort2_A1','PROJ1_Cohort2_A2', 'PROJ1_Cohort2_B1', 'PROJ1_Cohort2_B2', 'PROJ1_Cohort2_C1', 'PROJ1_Cohort2_C2',
              'PROJ1_Cohort3_A1', 'P1_Cohort3_A2')

# Array type, used : EPIC or 450K
artype <- '450K'
# Parameters to exclude CpGs
exclude <- c( 'MASK_sub30_copy', 'MASK_extBase', 'MASK_mapping', 'MASK_typeINextBaseSwitch', 'control_probes', 'Unrel_450_EPIC_blood', 'Sex')

# Ethnic group
ethnic <- 'EUR'

N <- c(100, 100, 166, 166, 166, 166, 166, 166, 240, 240 )
n <- c(NA)


# Venn diagrams
venn_diagrams <- list(
   c("Cohort1_A1", "PROJ1_Cohort2_A1", "PROJ1_Cohort2_B1", "PROJ1_Cohort2_C1", "PROJ1_Cohort3_A1" ),
   c("Cohort1_A2", "PROJ1_Cohort2_A2", "PROJ1_Cohort2_B2", "PROJ1_Cohort2_C2", "P1_Cohort3_A2" )
)

########## ----------  END VARIABLES DEFINED BY USER  ----------  ##########



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

   # Remove rows with NA from data
   cohort <- clean_NA_from_data(cohort)

   # Descriptives - Before CpGs deletion #
   descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, before = TRUE)

   # Remove duplicates
   # cohort <- remove_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i], '/',prefixes[i],'_descriptives_duplic.txt'), paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )

   test_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )

   # Exclude CpGs not meet conditions
   cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = ethnic, filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName )

   # Descriptives - After CpGs deletion #
   descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, before = FALSE )

   # Adjust data by Bonferroni and FDR
   cohort <- adjust_data(cohort, "P_VAL", bn=TRUE, fdr=TRUE, fResumeName  )

   # Write QC complete data to external file
   write_QCData(cohort, paste0(results_folder, '/',prefixes[i], '/',prefixes[i]))


   ## Visualization - Plots
   rasterpdf::raster_pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'), res = 300)

   # Distribution plot
   plot_distribution(cohort$SE, main = paste('Standard Errors of', prefixes[i]), xlab = 'SE')
   plot_distribution(cohort$P_VAL, main = paste('p-values of', prefixes[i]), xlab = 'p-value')

   # QQ plot.
   qqman::qq(cohort$P_VAL, main = sprintf('QQ plot of %s (lambda = %f)', prefixes[i], lambda=get_lambda(cohort,"P_VAL")))

   # Volcano plot.
   plot_volcano(cohort, "BETA", "P_VAL", main =paste('Volcano plot of', prefixes[i]) )

   dev.off()

   # Add mandatory data for precisionplot
   medianSE[i] <-  median(cohort$SE)
   value_N[i] <- N[i]
   cohort_label[i] <-  prefixes[i]

   # if n is defined for dichotomic condition :
   if(length(n) == length(N))  value_n[i] <- n[i]

   # Store data for Beta Box-Plot
   if( i == 1)
      betas.data <- list()
   betas.data[[prefixes[i]]] <- cohort[,"BETA"]

}

# Data for Precision Plot
precplot.data <- as.data.frame(cbind( SE = medianSE, invSE = (1/medianSE), N = value_N, sqrt_N = sqrt(N), cohort = cohort_label ))

if(length(n) == length(N))
   precplot.data.n <- as.data.frame(cbind( SE = medianSE, invSE = (1/medianSE), N = value_n, sqrt_N = sqrt(n), cohort = cohort_label ))

# BoxPlot with Betas in all Models and cohorts
plot_betas_boxplot(betas.data, paste(results_folder, 'BETAS_BoxPlot.pdf', sep="/"))

##  Post model analysis  ##

if ( length(files) > 1)
{
   # Precision_Plot(N)
   plot_precisionp(precplot.data, paste(results_folder,  "precision_SE_N.png",sep='/'), main = "Precision Plot - 1/median(SE) vs sqrt(N)")

   # Precision_Plot(n)
   if(length(n) == length(N))
      plot_precisionp(precplot.data.n, paste(results_folder,  "precision_SE_N.png", sep='/'), main = "Subgroup Precision Plot -  1/median(SE) vs sqrt(n)")

   # Venn_Diagrams()
   for (i in 1:length(venn_diagrams))
      plot_venndiagram(venn_diagrams[[i]], qcpath = results_folder, plotpath = results_folder, bn='padj.bonf', fdr='padj.fdr')

}
