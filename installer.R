# Install requirede libraries
devtools::install_github("ilarischeinin/rasterpdf")
devtools::install_github("guido-s/meta")
devtools::install_version('tidyverse', version = '1.3.0')
devtools::install_version('dplyr', version = '1.0.2')
devtools::install_version('stringr', version = '1.4.0')
devtools::install_version('ggplot2', version = '3.3.2')
devtools::install_version('VennDiagram', version = '1.6.20')
devtools::install_version('RColorBrewer', version = '1.1-2')
devtools::install_version('reshape2', version = '1.4.4')
devtools::install_version('ggsignif', version = '0.6.0')
install.packages("tools")

# Install brgeEnrich package from github-brge repositorie for enrichment :

# Install required packages for brgeEnrich
devtools::source_url("https://raw.githubusercontent.com/isglobal-brge/BrgeEnrich/HEAD/installer.R")
devtools::install_github("isglobal-brge/brgeEnrich@HEAD")

#..# devtools::install_version('tools', version = '3.6.3', repos = 'https://cran.us.r-project.org' )
devtools::install_version('BiocManager', version = '1.30.10')

# Bioconductor
BiocManager::install( "missMethyl", ask = FALSE )
BiocManager::install( "org.Hs.eg.db", ask = FALSE )
BiocManager::install( "GenomicRanges", ask = FALSE )
BiocManager::install( "rtracklayer", ask = FALSE )
BiocManager::install( "IlluminaHumanMethylation450kanno.ilmn12.hg19", ask = FALSE )
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", ask = FALSE)
