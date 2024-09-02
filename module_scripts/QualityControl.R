## ################################################## ##
##  Quality Control Script to use with EASIER package ##
##                                                    ##
##  script version; 0.1.2.27                          ##
## ################################################## ##


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


# load package
library(EASIER)

########## ----------  VARIABLES DEFINED BY USER  ----------  ##########

# (this data is only an example)

# Set working directory to metaanalysis folder
setwd("<path to metaanalysis folder>/metaanalysis")

# Files used in QC, needed in meta-analysis to plot ForestPlot
files <- c('data/PACE_GENR_20201006_C.txt',
           'data/PACE_GENR_20201006_C.txt',
           'data/Cohort1_Model1_20170713.txt',
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
prefixes <- c('PACE_GENR_C_EUR','PACE_GENR_C_GMAF1p',
              'Cohort1_A1', 'Cohort1_A2',
              'PROJ1_Cohort2_A1','PROJ1_Cohort2_A2', 'PROJ1_Cohort2_B1', 'PROJ1_Cohort2_B2', 'PROJ1_Cohort2_C1', 'PROJ1_Cohort2_C2',
              'PROJ1_Cohort3_A1', 'P1_Cohort3_A2')


# Exclude - MASK snp5
ethnic <- c('EUR','GMAF1p', 'EUR', 'SAS', 'EUR', 'EAS', 'EUR', 'SAS', 'EUR', 'EUR', 'EUR', 'EAS')

# Array type, used : EPIC or 450K
artype <- c('450K', '450K', 'EPIC', '450K', 'EPIC', '450K', '450K', '450K', 'EPIC', '450K', '450K', 'EPIC')

# Parameters to exclude CpGs
exclude <- c( 'control_probes',
              'noncpg_probes',
              'Sex',
              'MASK_mapping',
              'MASK_sub30_copy',
              'MASK_extBase',
              'MASK_typeINextBaseSwitch',
              'MASK_snp5_ethnic',
              'Unrel_450_EPIC_blood')



N <- c(100, 100, 100, 100, 166, 166, 166, 166, 166, 166, 240, 240 )
n <- c(NA)

# Minimum sample representation percentage required for CpGs
# Filter minimum percentage of missings for each CpG in cohort
# We need to define two parameters,
#  - colname_NforProbe:
#        Column name with Number of individuals per probe, this variable only needs
#           to be defined if you want to filter CpGs with low representation.
#         If defined value in colname_NforProbe not exists, no filter will be applied
#  - pcMissingSamples :
#        Máximum percent of missing samples allowed,

colname_NforProbe <- 'N_for_probe'
pcMissingSamples <- 0.9

# Venn diagrams
venn_diagrams <- list(
   c("Cohort1_A1", "PROJ1_Cohort2_A1", "PROJ1_Cohort2_B1", "PROJ1_Cohort2_C1", "PROJ1_Cohort3_A1" ),
   c("Cohort1_A1", "PROJ1_Cohort2_A1"),
   c("Cohort1_A2", "PROJ1_Cohort2_A2", "PROJ1_Cohort2_B2", "PROJ1_Cohort2_C2", "P1_Cohort3_A2" )
)

########## ----------  END VARIABLES DEFINED BY USER  ----------  ##########



## ###################### ##
##  QC - Quality Control  ##
## ###################### ##

# Variable declaration to perform precision plot
medianSE <- numeric(length(files))
value_N <- numeric(length(files))

if(length(n) == length(N))
   value_n <- numeric(length(files))

cohort_label <- character(length(files))

# Prepare output folder for results (create if not exists)
if(!dir.exists(file.path(getwd(), results_folder )))
   suppressWarnings(dir.create(file.path(getwd(), results_folder)))

# Remove tmp files
if( file.exists(paste0(results_folder,"/tmp_pretQC.txt")) )
    file.remove(paste0(results_folder,"/tmp_pretQC.txt"))
if( file.exists(paste0(results_folder,"/tmp_pretQC.txt")) )
    file.remove(paste0(results_folder,"/tmp_postQC.txt")) 
if( file.exists(paste0(results_folder,"/tmp_pretQC.txt")) )
    file.remove(paste0(results_folder,"/tmp_postQCAdj.txt")) 

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
   cohort <- descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], artype = artype[i], before = TRUE)

   # Remove duplicates
   # cohort <- remove_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i], '/',prefixes[i],'_descriptives_duplic.txt'), paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )

   test_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )

   # Remove cpGs with low representation
   # first, we test if colname_NforProbe and pcMissingSampes are defined
   if( !exists("colname_NforProbe") ) { colname_NforProbe <- NULL }
   if( !exists("pcMissingSamples") ) { pcMissingSamples <- NULL }

   cohort <- filterLowRepresentedCpGsinCohort(cohort, colname_NforProbe, pcMissingSamples, N[i], fileresume = fResumeName )

   # Exclude CpGs not meet conditions
    if("MASK_snp5_ethnic" %in% exclude ){
       cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = ethnic[i], filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
    } else {
       #..# if( !is.null(exclude) && exclude!='') {
         cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = "", filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
         #..# }
    }

   # Descriptives - After CpGs deletion #
   descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], artype = artype[i], before = FALSE )

   # Adjust data by Bonferroni and FDR
   cohort <- adjust_data(cohort, "P_VAL", bn=TRUE, fdr=TRUE, fResumeName, N[i]  )

   # Write QC complete data to external file
   write_QCData(cohort, paste0(results_folder, '/',prefixes[i], '/',prefixes[i]))

   ## Visualization - Plots
   #. Problems in some workstations and servers.# rasterpdf::raster_pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'), res = 300)
   #..# Problems in some cases --> Get png plots : pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'))

   # Distribution plot
   png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_SE_plot.png'), type="cairo")
      plot_distribution(cohort$SE, main = paste('Standard Errors of', prefixes[i]), xlab = 'SE')
   dev.off()
   png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_pvals_plot.png'), type="cairo")
      plot_distribution(cohort$P_VAL, main = paste('p-values of', prefixes[i]), xlab = 'p-value')
   dev.off()

   # QQ plot
   png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_QQ_plot.png'), type="cairo")
      qqman::qq(cohort$P_VAL, main = sprintf('QQ plot of %s (lambda = %f)', prefixes[i], lambda=get_lambda(cohort,"P_VAL")))
   dev.off()

   # Volcano plot.
   png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_Volcano_plot.png'), type="cairo")
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

# Create QC Summary
if( file.exists(paste0(results_folder,"/tmp_pretQC.txt")) & file.exists(paste0(results_folder,"/tmp_postQC.txt")) & file.exists(paste0(results_folder,"/tmp_postQCAdj.txt") )) {
   preQC <- read.table (file = paste0(results_folder,"/tmp_pretQC.txt"), header = TRUE, sep = "\t")
   postQC <- read.table (file = paste0(results_folder,"/tmp_postQC.txt"), header = TRUE, sep = "\t")
   postQCAdj <- read.csv(file = paste0(results_folder,"/tmp_postQCAdj.txt"), header = TRUE, sep = "\t")

   write.table( cbind(prefixes,ethnic, postQC[,1:2], preQC, postQC[,3:length(postQC)], postQCAdj), file = paste0(results_folder,"/Summary_QCs.txt" ),
                row.names = FALSE, col.names = TRUE, sep = "\t")

   do.call(file.remove, list(list.files(results_folder, full.names = TRUE, pattern = "tmp_*")))
}


# Data for Precision Plot
precplot.data <- cbind.data.frame( SE = medianSE, invSE = (1/medianSE), N = value_N, sqrt_N = sqrt(N), cohort = cohort_label )
cols.numeric <- c("SE","invSE", "N", "sqrt_N")
precplot.data[cols.numeric] <- sapply(precplot.data[cols.numeric],as.numeric)

if(length(n) == length(N)){
   precplot.data.n <- cbind.data.frame( SE = medianSE, invSE = (1/medianSE), N = value_n, sqrt_N = sqrt(n), cohort = cohort_label )
   precplot.data.n[cols.numeric] <- sapply(precplot.data.n[cols.numeric],as.numeric)
}

# BoxPlot with Betas in all Models and cohorts
plot_betas_boxplot(betas.data, paste(results_folder, 'BETAS_BoxPlot.png', sep="/"))

##  Post model analysis  ##

if ( length(files) > 1)
{
   # Precision_Plot(N)
   plot_precisionp(precplot.data, paste(results_folder,  "precision_SE_N.png",sep='/'), main = "Precision Plot - 1/median(SE) vs sqrt(n)")

   # Precision_Plot(n)
   if(length(n) == length(N))
      plot_precisionp(precplot.data.n, paste(results_folder,  "precision_SE_N.png", sep='/'), main = "Subgroup Precision Plot -  1/median(SE) vs sqrt(n)")

   # Venn_Diagrams()
   for (i in 1:length(venn_diagrams))
      plot_venndiagram(venn_diagrams[[i]], qcpath = results_folder, plotpath = results_folder, bn='padj.bonf', fdr='padj.fdr')

}
