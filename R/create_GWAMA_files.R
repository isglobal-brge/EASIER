#' Creates GWAMA formatted files
#'
#' Creates files formatted with GWAMA specifications to perform a meta-analysis with GWAMA software
#'
#' @param qcpath string. Path with Quality Control results (_QCData.txt files)
#' @param filename string. File name to formatt to use with GWAMA
#' @param gwmodelpath string. Path to store gwama meta-analisis files, all gwama files related to the same meta-analysis must be under the same path
#' @param N number. Number of samples in model
#' @param llowCpGs Optiona, numeric array with all CpGs with low representation in models.
#'
#' @return files with model data formatted for GWAMA
#'
#'
#' @export
create_GWAMA_files <- function(qcpath, file, gwmodelpath, N, llowCpGs=NULL)
{

   cohort <- read.table( paste0(qcpath,"/",file,"_QCData.txt")  , header = TRUE, as.is = TRUE, dec='.')
   cohort <- tibble::as_tibble(cohort)

   # If we have lowCpGs list, we filter this CpGs in the new file
   if( !is.null(llowCpGs) ){
      # Rmove CpGs with exclusion parameters
      cohort <- cohort[ !cohort$probeID %in% llowCpGs,]
   }

   # Prepare data for GWAMA
   data <- cohort %>% dplyr::select("MARKERNAME" = probeID, "CHR"= CpG_chrm, "POS"= CpG_beg, "BETA"= BETA, "SE"=SE)
   data <- data %>% tibble::add_column( STRAND = rep('+',dim(cohort)[1]), .after = "MARKERNAME")%>%
      tibble::add_column(EA = rep('C',dim(cohort)[1]), .after = "POS")%>%
      tibble::add_column(NEA = rep('G',dim(cohort)[1]), .after = "EA")%>%
      tibble::add_column(N = rep(N,dim(cohort)[1]), .after = "SE" )

   # To prevent changes from '.' separator in double to ','
   data %>%
      dplyr::mutate_at(vars(SE,BETA), character())
   data$BETA <- gsub("\\,", ".", data$BETA)
   data$SE <- gsub("\\,", ".", data$SE)

   qc.fname.gwamadata <- paste0(gwmodelpath, "/",file , '_gwama.txt')

   # Write file with GWAMA data
   suppressWarnings( write.table( data, qc.fname.gwamadata,
                                  col.names = TRUE, row.names = FALSE, sep = '\t', append = FALSE, quote = FALSE))


   # Add current file in gwama.in (gwama configuration file) to be processed by GWAMA
   qc.fname.gwamaconf <- paste(gwmodelpath,"gwama.in", sep = "/")
   qc.fname.gwamadata <- strsplit(qc.fname.gwamadata,"/")[[1]][length(strsplit(qc.fname.gwamadata,"/")[[1]])]
   qc.fname.gwamadata <- paste(getwd(), gwmodelpath, qc.fname.gwamadata,sep="/")

   if( i == 1 ) {
      write(sprintf('%s', qc.fname.gwamadata), file = qc.fname.gwamaconf)
   } else  {
      write(sprintf('%s', qc.fname.gwamadata), file = qc.fname.gwamaconf, append = TRUE) }


}
