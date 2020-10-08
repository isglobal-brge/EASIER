#' Forest plot
#'
#' Forest plot
#'
#' @param datas complete data obtained from meta-analysis after QC without annotations
#' @param files_meta string vector with files used in the meta-analysis
#' @param islowCpg string  string indicating if the execution has been carried out with all the CpGs or with the filtered CpGs, possible values : 'Normal', 'lowcpgs'
#' @param gwana_dir string with gwama input data path
#' @param metaname string with meta-analysis name
#'
#' @return distribution plot
#'
#' @export
plot_ForestPlot <- function( datas, files_meta, islowCpg, gwana_dir, metaname )
{
   # We need the data from meta-analysis model and individual data from cohorts in the model

   ## FER-HO AMB UN APPLY

   type = c('Fixed','Random')

   names(files_meta) <- files_meta
   cohorts_names <- split(unname(files_meta),names(files_meta))

   # Read data from pre-meta-analysis cohorts and sort by p-value
   cohorts <- lapply(files_meta, function(cohortn) {
      #..# cohortfile <- paste0(gwana_dir,"/", cohortn,"_gwama.txt")
      cohortfile <- files[which(prefixes==cohortn)]
      print(paste0("Llegim : ",cohortfile))
      cohortdata <- read.table(cohortfile, header = TRUE, stringsAsFactors=FALSE )
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

      data <- data[order(data$p.value),]

      # We only get the 30 first CpGs with lowest p-value
      cpgs <- as.character(data[which(data$p.value<0.05),'rs_number'][1:30])

      tt <- data.frame(do.call(rbind, lapply(cpgs, function(cp) {
         tt <- sapply(cohorts, function(ch) {
            x <- which(ch$probeID == cp)
            if( length( x ) == 0) { NA }
            else { x }
         })
         names(tt) <- names(cohorts)
         tt
      })))
      rownames(tt) <- cpgs

      bb <- lapply(rownames(tt), function(cpg) {
         data.frame(t(sapply(names(cohorts), function(ch) {
            cohorts[[ch]][tt[cpg, ch][[1]], 2:3]
         })))
      })

      names(bb) <- rownames(tt)

      if(islowCpg == 'lowcpgs')
         path <- paste0(file.path(getwd()),"/", outputfolder,"/",metaname,"_Filtr/ForestPlots")
      else
         path <- paste(file.path(getwd()), outputfolder,metaname,"ForestPlots", sep="/")

      # Before get plots test if dir exists and create it
      if(!dir.exists(path))
         suppressWarnings(dir.create(path))

      mt <- lapply(names(bb), function(cpg) {
         message(cpg)
         df <- bb[[cpg]]
         mtg <- meta::metagen(unlist(df[,"BETA"]), unlist(df[,"SE"]), sm="MD", studlab=rownames(df), comb.random=TRUE, comb.fixed = TRUE)

         #pdf(paste0( path, "/FP_", cpg,"_",type[ts] ,".pdf"))
         #rasterpdf::
         rasterpdf::raster_pdf(paste0( path, "/FP_", cpg,"_",type[ts] ,".pdf"), res = 600)
         par(mar = c(0, 0, 0, 0))
         meta::forest(mtg, leftcols=c("studlab"), leftlabs=c("Cohort"), rightcols=c("effect", "ci","pval","w.fixed","w.random"), fontsize=7, digits=3, print.pval=TRUE, addrow.overall=T,
                col.fixed="red", col.random="blue",print.tau2 = FALSE, smlab = "", col.diamond.fixed="red", col.diamond.random = "blue", overall= T, test.overall=T,
                fs.test.overall=7, fs.hetstat=5, fs.axis=5)
         dev.off()
      })

   }

}