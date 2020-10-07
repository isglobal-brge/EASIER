if (!require(rasterpdf, quietly = TRUE)) install.packages('rasterpdf')
if (!require(ggplot2, quietly = TRUE)) install.packages('ggplot2')
if (!require(VennDiagram, quietly = TRUE)) install.packages('VennDiagram')
if (!require(RColorBrewer, quietly = TRUE)) install.packages('RColorBrewer')
library(methyTools)

#..# setwd("~/Library/Mobile Documents/com~apple~CloudDocs/PROJECTES/Treballant/methyTools")
setwd("/Users/mailos/Library/Mobile Documents/com~apple~CloudDocs/ISGlobal/Projectes/EWAS/2_Llibreria/proves")

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

#

## ############### ##
##  Meta-Analysis  ##
## ############### ##





