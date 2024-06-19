# Install requirede libraries
if (!require('qqman')){ install.packages("qqman") }
if (!require('rasterpdf')){ devtools::install_github("ilarischeinin/rasterpdf") }
if (!require('meta')){ devtools::install_github("guido-s/meta") }
if (!require('tidyverse')){ install.packages('tidyverse') }
if (!require('dplyr')){ install.packages('dplyr') }
if (!require('stringr')){ install.packages('stringr') }
if (!require('ggplot2')){ install.packages('ggplot2') }
if (!require('VennDiagram')){ install.packages('VennDiagram') }
if (!require('RColorBrewer')){ install.packages('RColorBrewer') }
if (!require('reshape2')){ install.packages('reshape2') }
if (!require('ggsignif')){ install.packages('ggsignif') }
if (!require('readtext')){ install.packages('readtext') }
if (!require('tools')){ install.packages("tools") }

# Install brgeEnrich package from github-brge repositorie for enrichment :

# Install required packages for brgeEnrich
devtools::source_url("https://raw.githubusercontent.com/isglobal-brge/BrgeEnrich/HEAD/installer.R")
devtools::install_github("isglobal-brge/brgeEnrich@HEAD")

if (!require("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager") 

# Bioconductor
if (!require('missMethyl')){ BiocManager::install( "missMethyl", ask = FALSE ) }
if (!require('org.Hs.eg.db')){ BiocManager::install( "org.Hs.eg.db", ask = FALSE ) }
if (!require('GenomicRanges')){ BiocManager::install( "GenomicRanges", ask = FALSE ) }
if (!require('rtracklayer')){ BiocManager::install( "rtracklayer", ask = FALSE ) }
if (!require('IlluminaHumanMethylation450kanno.ilmn12.hg19')){ BiocManager::install( "IlluminaHumanMethylation450kanno.ilmn12.hg19", ask = FALSE ) }
if (!require('IlluminaHumanMethylationEPICanno.ilm10b4.hg19')){ BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", ask = FALSE) }
