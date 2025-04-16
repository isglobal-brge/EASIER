#' 450K Illumina filters
#'
#' A dataset containing the 450K Illumina methylation array probe filtering conditions
#'
#' @docType data
#'
#' @usage data(filter_450K)
#'
#' @format A data frame with 485577 rows and 50 variables
#' @references Chen et al. 2013
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
#'   \item{Unrel_450_EPIC_blood }{Unreliable probes discordant between 450K and EPIC for blood (TRUE/FALSE). }
#'   \item{Unrel_450_EPIC_pla }{Unreliable probes discordant between 450K and EPIC for placenta (TRUE/FALSE). }
#'   \item{Unrel_450_EPIC_pla_restrict }{Unreliable probes discordant between 450K and EPIC for placenta, more restrictive (TRUE/FALSE). }
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
#' @format A data frame with 865918 rows and 50 variables
#' @references Chen et al. 2013
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
#'   \item{Unrel_450_EPIC_blood }{Unreliable probes discordant between 450K and EPIC for blood (TRUE/FALSE). }
#'   \item{Unrel_450_EPIC_pla }{Unreliable probes discordant between 450K and EPIC for placenta (TRUE/FALSE). }
#'   \item{Unrel_450_EPIC_pla_restrict }{Unreliable probes discordant between 450K and EPIC for placenta, more restrictive (TRUE/FALSE). }
#'   \item{CpG_chrm}{ Chromosome }
#' }
"filter_EPIC"

#' Molecular Signatures Database (MSigDB v 7.1) - C2 curated gene sets
#'
#' A dataset containing Molecular Signatures Database in his current version 7.1.
#'
#' @docType data
#'
#' @usage data(human_c2_v7p1)
#'
#' @format A data frame with 5529 variables
#' @references \url{http://bioinf.wehi.edu.au/MSigDB/v7.1/}
"Hs.c2"


#' Chromatine states data
#'
#' Chromatine states data
#'
#' @docType data
#'
#' @usage data(crom15)
#' @format A data frame with 482421 CpGs and 21 variables related to Chromatine states
#' @references BRGE
#' \describe{
#'		\item{X}{X}
#'		\item{HT12v4.ArrayAddress}{HT12v4.ArrayAddress}
#'		\item{Gene}{Gene}
#'		\item{Chr}{Chr}
#'		\item{ChrStart}{ChrStart}
#'		\item{ChrEnd}{ChrEnd}
#'		\item{Probe}{Probe}
#'		\item{TssA}{TssA}
#'		\item{TssAFlnk}{TssAFlnk}
#'		\item{TxFlnk}{TxFlnk}
#'		\item{TxWk}{TxWk}
#'		\item{Tx}{Tx}
#'		\item{EnhG}{EnhG}
#'		\item{Enh}{Enh}
#'		\item{ZNF.Rpts}{ZNF.Rpts}
#'		\item{Het}{Het}
#'		\item{TssBiv}{TssBiv}
#'		\item{BivFlnk}{BivFlnk}
#'		\item{EnhBiv}{EnhBiv}
#'		\item{ReprPC}{ReprPC}
#'		\item{ReprPCWk}{ReprPCWk}
#'		\item{Quies}{Quies}
#' }
"crom15"

#' dhs data
#'
#' dhs data
#'
#' @docType data
#'
#' @usage data(dhs)
#' @format A data frame with 482421 CpGs and 19 variables
"dhs"

#' Fetal Placenta 15 Stats
#'
#' Fetal placentat 15 states (E091)
#'
#' @docType data
#'
#' @usage data(FP_15_E091)
#' @format Genomic Ranges
"FP_15_E091"

#' Fetal Placenta 18 Stats
#'
#' Fetal placentat 18 states (E091)
#'
#' @docType data
#'
#' @usage data(FP_18_E091)
#'
#' @format Genomic Ranges
"FP_18_E091"


#' PMD placenta
#'
#' PMD placenta
#'
#' @docType data
#'
#' @usage data(PMD_placenta)
#' @references Schroeder, D. I. et al. The human placenta methylome. Proc. Natl. Acad. Sci. 110, 6037–6042 (2013)
#' @format A data frame with 3 variables
"PMD_placenta"

#' Imprinting Regions Placenta
#'
#' Imprinting regions for placenta
#'
#' @docType data
#'
#' @usage data(IR_Placenta)
#' @format A data frame with  regions and 5 variables
#' \describe{
#'   \item{Chr_DMR}{Chromosome}
#'   \item{Start_DMR}{Start position}
#'   \item{End_DMR}{End position}
#'   \item{DMR_Classification}{DMR classification, possible values are "Candidate mDMR", "Known placenta-specific mDMR", "Candidate pDMR", "Known mDMR" and "Known pDMR".}
#'   \item{Known_gDMR }{Known DMR associated gene }
#' }
"IR_Placenta"

#' Imprinting Regions Carrera Epic V1
#'
#' Imprinting regions for Epic Version 1
#'
#' @docType data
#'
#' @usage data(IRCarrera_Ev1)
#' @format A data frame with 156 regions
#' \describe{
#'   \item{ICR_id}{ICR_id}
#'   \item{CpG_chr}{Chromosome}
#'   \item{CpG_start}{Start position}
#'   \item{CpG_stop}{End position}
#'   \item{ICR_chr}{Chromosome}
#'   \item{ICR_start}{Start position}
#'   \item{ICR_stop}{End position}
#'   \item{CpG_Probe}{DMR classification, possible values are "Candidate mDMR", "Known placenta-specific mDMR", "Candidate pDMR", "Known mDMR" and "Known pDMR".}
#'   \item{ID}{Original ICR ID}
#'   \item{Genom_Coord}{Genomic Coordinates}
#'   \item{Meth_Origin}{Parental Origin of Methylation}
#'   \item{Near_Trasncript}{Nearest Transcript}
#'   \item{Dist_Near_Trasncript}{Distance to Nearest Transcript}
#' }
"IRCarrera_Ev1"


#' Imprinting Regions Carrera Epic V2
#'
#' Imprinting regions for Epic Version 2
#'
#' @docType data
#'
#' @usage data(IRCarrera_Ev2)
#' @format A data frame with 548 regions
#' \describe{
#'   \item{ICR_id}{ICR_id}
#'   \item{CpG_chr}{Chromosome}
#'   \item{CpG_start}{Start position}
#'   \item{CpG_stop}{End position}
#'   \item{ICR_chr}{Chromosome}
#'   \item{ICR_start}{Start position}
#'   \item{ICR_stop}{End position}
#'   \item{CpG_Probe}{DMR classification, possible values are "Candidate mDMR", "Known placenta-specific mDMR", "Candidate pDMR", "Known mDMR" and "Known pDMR".}
#'   \item{ID}{Original ICR ID}
#'   \item{Genom_Coord}{Genomic Coordinates}
#'   \item{Meth_Origin}{Parental Origin of Methylation}
#'   \item{Near_Trasncript}{Nearest Transcript}
#'   \item{Dist_Near_Trasncript}{Distance to Nearest Transcript}
#' }
"IRCarrera_Ev2"

#' eQTM filter
#'
#' Data for eQTM enrichment
#'
#' @docType data
#'
#' @usage data(eQTM_autosome_unadj_filtered.cells)
#' @format A data frame with 104716 CpGs and 4 variables
#' \describe{
#'   \item{CpG}{CpG id}
#'   \item{p.value}{p-value}
#'   \item{sigPair}{is a significant Pair}
#'   \item{TC_gene}{Gene related to CpG}
#' }
"eQTM"
