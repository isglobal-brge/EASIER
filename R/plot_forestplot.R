#' Forest plot
#'
#' Forest plot
#'
#' @param datas complete data obtained from meta-analysis after QC without annotations
#' @param files_meta string vector with files used in the meta-analysis
#' @param islowCpg string  string indicating if the execution has been carried out with all the CpGs or with the filtered CpGs, possible values : 'Normal', 'lowcpgs'
#' @param gwana_dir string with gwama input data path
#' @param metaname string with meta-analysis name
#' @param files string vector with files used in QC because we need data from this files to perform ForestPlots
#' @param outputgwama string with gwama output path
#' @param nsignificatives number, lowest p-values to show in forestplot, by default, lowestp=30
#' @param criteria character string with criteria to be used to select the first significative CpGs,
#'  by default the criteria is p.value (without any adjustment). Possible values are "p.value", "FDR" or "Bonferroni", in that case function plots all significative CpGs.
#' @param selectedCpGs character vector with the names of the CpGs to be plotted
#'
#' @return distribution plot
#' @importFrom data.table like
#' @export
plot_ForestPlot <- function( datas, files_meta, islowCpg, gwana_dir, metaname, files, outputgwama, nsignificatives = 30, criteria = "p.value" , selectedCpGs = NULL )
{
   # We need the data from meta-analysis model and individual data from cohorts in the model

   if(criteria )


   if(is.null(selectedCpGs)) {
      available.criteria <- c( "p.value", "FDR", "Bonferroni")
      method <- match(criteria, available.criteria)

      if (any(is.na(method))) {
         stop("You wrote the name of an unavailable criteria Available criteria are:
               p.value, FDR and Bonferroni")
      }
   }

   type = c('Fixed','Random')

   names(files_meta) <- files_meta
   cohorts_names <- split(unname(files_meta),names(files_meta))

   # Read data from pre-meta-analysis cohorts and sort by p-value
   cohorts <- lapply(files_meta, function(cohortn) {
      #..# cohortfile <- paste0(gwana_dir,"/", cohortn,"_gwama.txt")
      cohortfile <- files[which(prefixes==cohortn)]
      print(paste0("Reading : ",cohortfile))
      cohortdata <- read.table(cohortfile, header = TRUE, stringsAsFactors=FALSE )
      if(ncol(cohortdata)==1) {
         cohortdata <- read.table(cohortfile, header = TRUE, stringsAsFactors=FALSE, sep = ',' )
      }
      if(ncol(cohortdata)==1) {
         stop(paste0("Error reading input file ", cohortfile))
      }
      colnames(cohortdata) <- toupper(colnames(cohortdata))
      cohortdata <- cohortdata[order(cohortdata$P_VAL),]
   })
   names(cohorts) <- files_meta
   ts <- 0
   # Get Forest plot for Fixed and Random effects
   for ( d in 1:length(datas) )    # ==> TOT AIXÒ S'HA DE VECTORITZAR
   {
      data <- datas[[d]]
      ts <- ts + 1
      # Sort data before get rs_numbers --> forestplot only for first 30 CpGs

      if(is.null(selectedCpGs)){

         if(criteria %in% c('p.value', "FDR")){
            data <- data[order(data[,criteria]),]
            # We only get the nsignificatives first CpGs with lowest p-value
            # cpgs <- as.character(data[which(data$p.value<0.05),'rs_number'][1:nsignificatives])
            nsignificatives <- length(data[which(data[,criteria]<0.05),'rs_number'])
            if( nsignificatives > nsignificatives ){
               cpgs <- as.character(data[which(data[,criteria]<0.05),'rs_number'][1:nsignificatives])
            } else {
               if(nsignificatives > 0){
                  cpgs <- as.character(data[which(data[,criteria]<0.05),'rs_number'][1:nsignificatives])
               } else {
                  cpgs <- NULL
               }
            }

         } else {
            cpgs <- as.character(data[which(data[,criteria]!="no"),'rs_number'])
         }

      } else {
         if(inherits(selectedCpGs, "character")){
            cpgs <- selectedCpGs
         } else {
            warning(paste0("selectedCpgs must be a character, please set selectedCpGs as character or NULL to plot first xxx signiifcative CpGs"))
            cpgs <- NULL
         }
      }

      if(!is.null(cpgs)) {

         tt <- data.frame(do.call(rbind, lapply(cpgs, function(cp) {
            tt <- sapply(cohorts, function(ch) {
               x <- which(ch[, which( toupper(colnames(ch)) == 'PROBEID' )]  == cp)
               if( length( x ) == 0) { NA }
               else { x }
            })
            names(tt) <- names(cohorts)
            tt
         })))
         rownames(tt) <- cpgs

         bb <- lapply(rownames(tt), function(cpg) {
            data.frame(t(sapply(names(cohorts), function(ch) {
               cohorts[[ch]][tt[cpg, ch][[1]], c("BETA","SE")]
            })))
         })

         names(bb) <- rownames(tt)

         outputfolder <- ifelse( substr(outputgwama, 1, 2) == './', substr(outputgwama,3,nchar(outputgwama)), outputgwama)

         if(islowCpg == 'lowcpgs') {
            path <- paste0(file.path(getwd()),"/", outputfolder,"/ForestPlots")
         } else {
            path <- paste(file.path(getwd()), outputfolder, "ForestPlots", sep="/")
         }

         # Before get plots test if dir exists and create it
         if(!dir.exists(path))
            suppressWarnings(dir.create(path, recursive = TRUE))

         mt <- lapply(names(bb), function(cpg) {
            message(cpg)
            dataf <- bb[[cpg]]
            if(packageVersion("meta") > "5.0.0") {
               mtg <- meta::metagen(unlist(dataf[,"BETA"]), unlist(dataf[,"SE"]), sm="MD", studlab=rownames(dataf), random=TRUE, fixed = TRUE)
            } else {
               mtg <- meta::metagen(unlist(dataf[,"BETA"]), unlist(dataf[,"SE"]), sm="MD", studlab=rownames(dataf), comb.random=TRUE, comb.fixed = TRUE)
            }


            # print(paste0("Output file : ",paste0( path, "/FP_", cpg,"_",type[ts] ,".pdf")))

            # rasterpdf::raster_pdf(paste0( path, "/FP_", cpg,"_",type[ts] ,".pdf"), res = 600)
            pdf(paste0( path, "/FP_", cpg,"_",type[ts] ,".pdf"))
            par(mar = c(0, 0, 0, 0))
            meta::forest(mtg, leftcols=c("studlab"), leftlabs=c("Cohort"), rightcols=c("effect", "ci","pval","w.fixed","w.random"), fontsize=7, digits=3, print.pval=TRUE, addrow.overall=T,
                         col.fixed="red", col.random="blue",print.tau2 = FALSE, smlab = "", col.diamond.fixed="red", col.diamond.random = "blue", overall= T, test.overall=T,
                         fs.test.overall=7, fs.hetstat=5, fs.axis=5, pooled.totals=TRUE)
            dev.off()
         })

      } else {
         warning(paste0("No significative CpGs found for ", criteria," criteria"))
      }



   }

}
