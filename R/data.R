#' 450K Illumina filters
#'
#' A dataset containing the 450K Illumina methylation array probe filtering conditions
#'
#' @docType data
#'
#' @usage data(filter_450K)
#'
#' @format A data frame with 485577 rows and 48 variables
#'
#' @references Chen et al. (2013)
#' \describe{
#'   \item{MASK.sub25.copy, MASK.sub30.copy, MASK.sub35.copy, MASK.sub40.copy}{indicate whether the 25bp, 30bp, 35bp and 40bp 3'-subsequence of the probe is non-unique (TRUE/FALSE)}
#'   \item{MASK.mapping}{"hether the probe is masked for mapping reason. Probes retained should have high quality (>=40 on 0-60 scale) consistent (with designed MAPINFO) mapping (for both in the case of type I) without INDELs (TRUE/FALSE). }
#'   \item{MASK.extBase }{Probes masked for extension base inconsistent with specified color channel (type-I) or CpG (type-II) based on mapping (TRUE/FALSE). }
#'   \item{MASK.typeINextBaseSwitch }{Whether the probe has a SNP in the extension base that causes a color channel switch from the official annotation (described as color-channel-switching, or CCS SNP in the reference). These probes should be processed differently than designed (by summing up both color channels instead of just the annotated color channel)(TRUE/FALSE).}
#'   \item{MASK.snp5.common }{Whether 5bp 3'-subsequence (including extension for typeII) overlap with any of the common SNPs from dbSNP (global MAF can be under 1%). }
#'   \item{MASK.snp5.GMAF1p }{Whether 5bp 3'-subsequence (including extension for typeII) overlap with any of the SNPs with global MAF >1% (TRUE/FALSE). }
#'   \item{MASK.snp5_<ethnicity> }{ One field for each possible ethnicyty : EUR SAS AMR GWD YRI TSI  IBS CHS PUR JPT  GIH CH_B STU ITU LWK KHV FIN ESN CEU PJL AC_B CLM CDX GBR BE_B PEL MSL  MXL ASW or  GMAF1p if population is very diverse. (ex:"MASK_snp5_EUR"). Whether 5bp 3'-subsequence (including extension for typeII) overlap with any of the SNPs for this particular ethnicity with a MAF >1% (TRUE/FALSE).}
#'   \item{MASK.general}{ Recommended general purpose masking merged from "MASK.sub30.copy", "MASK.mapping", "MASK.extBase", "MASK.typeINextBaseSwitch" and "MASK.snp5.GMAF1p" (TRUE/FALSE). }
#'   \item{probeType}{ cpg probes, non-cpg probes, and control probes are classified as "cg", "ch" or "rs" respectively in the variable named "probeType". }
#'   \item{Unreliable_450_EPIC }{Unreliable probes discordant between 450K and EPIC (TRUE/FALSE). }
#'   \item{CpG_chrm}{ Chromosome}
#' }
"filter_450K"


#' EPIC Illumina filters
#'
#' A dataset containing the EPIC Illumina methylation array probe filtering conditions
#'
#' @docType data
#'
#' @usage data(filter_EPIC)
#'
#' @format A data frame with 865918 rows and 48 variables
#' @references Chen et al. (2013)
#' \describe{
#'   \item{MASK.sub25.copy, MASK.sub30.copy, MASK.sub35.copy, MASK.sub40.copy}{indicate whether the 25bp, 30bp, 35bp and 40bp 3'-subsequence of the probe is non-unique (TRUE/FALSE)}
#'   \item{MASK.mapping}{"hether the probe is masked for mapping reason. Probes retained should have high quality (>=40 on 0-60 scale) consistent (with designed MAPINFO) mapping (for both in the case of type I) without INDELs (TRUE/FALSE). }
#'   \item{MASK.extBase }{Probes masked for extension base inconsistent with specified color channel (type-I) or CpG (type-II) based on mapping (TRUE/FALSE). }
#'   \item{MASK.typeINextBaseSwitch }{Whether the probe has a SNP in the extension base that causes a color channel switch from the official annotation (described as color-channel-switching, or CCS SNP in the reference). These probes should be processed differently than designed (by summing up both color channels instead of just the annotated color channel)(TRUE/FALSE).}
#'   \item{MASK.snp5.common }{Whether 5bp 3'-subsequence (including extension for typeII) overlap with any of the common SNPs from dbSNP (global MAF can be under 1%). }
#'   \item{MASK.snp5.GMAF1p }{Whether 5bp 3'-subsequence (including extension for typeII) overlap with any of the SNPs with global MAF >1% (TRUE/FALSE). }
#'   \item{MASK.snp5_<ethnicity> }{ One field for each possible ethnicyty : EUR SAS AMR GWD YRI TSI  IBS CHS PUR JPT  GIH CH_B STU ITU LWK KHV FIN ESN CEU PJL AC_B CLM CDX GBR BE_B PEL MSL  MXL ASW or  GMAF1p if population is very diverse. (ex:"MASK_snp5_EUR"). Whether 5bp 3'-subsequence (including extension for typeII) overlap with any of the SNPs for this particular ethnicity with a MAF >1% (TRUE/FALSE).}
#'   \item{MASK.general}{ Recommended general purpose masking merged from "MASK.sub30.copy", "MASK.mapping", "MASK.extBase", "MASK.typeINextBaseSwitch" and "MASK.snp5.GMAF1p" (TRUE/FALSE). }
#'   \item{probeType}{ cpg probes, non-cpg probes, and control probes are classified as "cg", "ch" or "rs" respectively in the variable named "probeType". }
#'   \item{Unreliable_450_EPIC }{Unreliable probes discordant between 450K and EPIC (TRUE/FALSE). }
#'   \item{CpG_chrm}{ Chromosome }
#' }
"filter_EPIC"

