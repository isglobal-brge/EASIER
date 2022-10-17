# Install requirede libraries
if (!require('qqman')){ install.packages("qqman") }
if (!require('rasterpdf')){ devtools::install_github("ilarischeinin/rasterpdf") }
if (!require('meta')){ devtools::install_github("guido-s/meta") }
if (!require('tidyverse')){ devtools::install_version('tidyverse', version = '1.3.0') }
if (!require('dplyr')){ devtools::install_version('dplyr', version = '1.0.2') }
if (!require('stringr')){ devtools::install_version('stringr', version = '1.4.0') }
if (!require('ggplot2')){ devtools::install_version('ggplot2', version = '3.3.2') }
if (!require('VennDiagram')){ devtools::install_version('VennDiagram', version = '1.6.20') }
if (!require('RColorBrewer')){ devtools::install_version('RColorBrewer', version = '1.1-2') }
if (!require('reshape2')){ devtools::install_version('reshape2', version = '1.4.4') }
if (!require('ggsignif')){ devtools::install_version('ggsignif', version = '0.6.0') }
if (!require('readtext')){ install.packages('readtext') }
if (!require('tools')){ install.packages("tools") }

# Install brgeEnrich package from github-brge repositorie for enrichment :

# Install required packages for brgeEnrich
devtools::source_url("https://raw.githubusercontent.com/isglobal-brge/BrgeEnrich/HEAD/installer.R")
devtools::install_github("isglobal-brge/brgeEnrich@HEAD")

#..# devtools::install_version('tools', version = '3.6.3', repos = 'https://cran.us.r-project.org' )
#..# devtools::install_version('BiocManager', version = '1.30.10')
if (!require("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager") 

# Bioconductor
if (!require('missMethyl')){ BiocManager::install( "missMethyl", ask = FALSE ) }
if (!require('org.Hs.eg.db')){ BiocManager::install( "org.Hs.eg.db", ask = FALSE ) }
if (!require('GenomicRanges')){ BiocManager::install( "GenomicRanges", ask = FALSE ) }
if (!require('rtracklayer')){ BiocManager::install( "rtracklayer", ask = FALSE ) }
if (!require('IlluminaHumanMethylation450kanno.ilmn12.hg19')){ BiocManager::install( "IlluminaHumanMethylation450kanno.ilmn12.hg19", ask = FALSE ) }
if (!require('IlluminaHumanMethylationEPICanno.ilm10b4.hg19')){ BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", ask = FALSE) }
