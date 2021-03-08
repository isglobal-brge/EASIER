#' Get descriptives after GWAMA
#'
#' Get descriptives after perform GWAMA meta-analysis
#'
#' @param resdir string. Path with results from GWAMA
#' @param analyzedata string vector. Route to analized data in order to get descriptives
#' @param modelfiles string vector with models used to perform meta-analysis
#' @param metaname string. Meta-analisis name. only used in titles
#' @param artype string, Illumina array type, 450K or EPIC
#'
#' @return File with descriptive and plots from GWAMA results. Results are stored in the same file as GWAMA results
#'
#'
#' @export
get_descriptives_postGWAMA <- function(resdir, analyzedata, modelfiles, metaname, artype, Ns)
{

   # Descriptives

   type <-  c('Fixed','Random')
   modeldata <- list()

   for( f in 1:length(analyzedata))
   {

      filedata <- paste0(analyzedata[f],".out")
      data <- read.table(filedata, header = TRUE, stringsAsFactors=FALSE )
      nCpG <- dim(data)[1]

      #..# qc.fname <- str_replace(analyzedata[f], type[f], "descriptives.txt")
      qc.fname <- paste0(analyzedata[f], "_descriptives.txt")

      print(paste0("Output file : ",qc.fname))

      write(sprintf('\t\t\t\t======================\n\t\t\t\t  Descriptive Resume\n\t\t\t\t======================\n'), file = qc.fname)
      write(sprintf('-----------------------------\n Model : %s\n-----------------------------\n',metaname), file = qc.fname, append = TRUE)
      write(sprintf('n models : %d ',length(unique(data$n_studies))), file = qc.fname, append = TRUE)
      write(sprintf('Models analysed :'), file = qc.fname, append = TRUE)
      write(sprintf('\t\t %s', unname(modelfiles)) , file = qc.fname, append = TRUE)
      write(sprintf('N Samples : %d ', sum(Ns) ), file = qc.fname, append = TRUE)
      write(sprintf('N CpGs : %d ',nCpG), file = qc.fname, append = TRUE)

      # Effects
      write(sprintf('\nEffects \n-------------\n'), file = qc.fname, append = TRUE)
      lst <- strsplit(data$effects, "") # Split effects in columns
      data <- cbind(data, data.frame(matrix(unlist(lst), nrow=length(lst), byrow=T)))
      effectpositions <- grep("X[0-9]",names(data)) # Posicions
      colnames(data)[grep("X[0-9]",names(data)) ] <- modelfiles   # effects colnames

      write(sprintf('# Effect from each file in model (relative freq.)\n'), file = qc.fname, append = TRUE)
      presence.pcent <- prop.table(do.call(cbind, lapply(data[effectpositions], table)), margin = 2)
      suppressWarnings(write.table(presence.pcent, qc.fname, sep="\t",row.names=TRUE, append = TRUE))

      # Descriptives
      descfields <- c(7,8,11,12,13,14,15,16)
      write(sprintf('\n# Descriptives for  Beta, SE, z, p-value, log10(p.value) ... , \n'), file = qc.fname, append = TRUE)
      desc <- prop.table(do.call(cbind, lapply(data[descfields], summary)), margin = 2)
      suppressWarnings(write.table( round(desc,5), qc.fname, sep="      \t",row.names=TRUE, append = TRUE))

      # Significative CpGs before and after Bonferroni Correction and FDR
      data$Bonferroni<-ifelse((data$p.value<0.05/dim(data)[1] ),"yes","no")   # Add Bonferroni correction addjustment
      data$FDR<-p.adjust(data$p.value, method="fdr")   # Add FDR adjustment

      write(sprintf('\n# Significative CpGs : \n'), file = qc.fname, append = TRUE)
      write(sprintf('\t# p-val<0.05 : %d\n', length(which(data$p.value<0.05)) ), file = qc.fname, append = TRUE)
      write(sprintf('\t# After Bonferroni Correction :  %d\n', length(which(data$Bonferroni == 'yes'))), file = qc.fname, append = TRUE)
      write(sprintf('\t# After FDR adjustmnet : %d\n', length(which(data$FDR<0.05)) ), file = qc.fname, append = TRUE)

      # Annotate CpGs
      data.annot <- annotate_CpGs(data[,"rs_number"], artype)

      # Merge data with annotation before write results to a file
      data.full <-  merge( data, data.annot, by.x = "rs_number", by.y ="Name")

      # Change "rs_number" column name by CpGId
      colnames(data.full)[1] <- "CpGId"

      # Remove allele columns
      data.full <- data.full[-grep("*allele",colnames(data.full))]

      # Write data with adjustments in other file, we don't want to overwrite gwama results
      qc.fname <- paste0(analyzedata[f], "_Modif.out")
      suppressWarnings(write.table(data.full, qc.fname, sep="\t",row.names=TRUE, append = FALSE))


      # QC Plots - Post metha-analysis

      # Custom graphical options.
      gg_colors <- function(n, alpha = NULL) {
         hues = seq(15, 375, length = n + 1)
         hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
      }




      # Distribution plot
      png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_SE_plot.png'))
      plot_distribution(cohort$SE, main = paste('Standard Errors of', prefixes[i]), xlab = 'SE')
      dev.off()
      png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_pvals_plot.png'))
      plot_distribution(cohort$P_VAL, main = paste('p-values of', prefixes[i]), xlab = 'p-value')
      dev.off()

      # QQ plot
      png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_QQ_plot.png'))
      qqman::qq(cohort$P_VAL, main = sprintf('QQ plot of %s (lambda = %f)', prefixes[i], lambda=get_lambda(cohort,"P_VAL")))
      dev.off()

      # Volcano plot.
      png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_Volcano_plot.png'))
      plot_volcano(cohort, "BETA", "P_VAL", main =paste('Volcano plot of', prefixes[i]) )
      dev.off()



      rasterpdf::raster_pdf(paste0(analyzedata[f], '_QCplots.pdf'), res = 600)

      # Plot distributions.
      plot.distr <- function(x, main, xlab) {
         h <- hist(x, breaks = 100, plot = FALSE)
         d <- density(x)
         plot(h, freq = FALSE, col = gg_colors(2)[2], border = 'white',
              main = main, xlab = xlab, ylim = c(0, max(d$y, h$density)))
         lines(d, col = gg_colors(2)[1], lwd = 2)
      }

      png(paste0(analyzedata[f], '_QC_distr_i2_plot.png'))
         plot.distr(na.omit(data$i2), main = paste('Heterogeneity (i2) histogram -', metaname,type[f]), xlab = 'i2')
      dev.off()

      png(paste0(analyzedata[f], '_QC_distr_SE_plot.png'))
         plot.distr(data$se, main = paste('Standard Errors -', metaname,type[f]), xlab = 'SE')
      dev.off()

      png(paste0(analyzedata[f], '_QC_distr_pvalue_plot.png'))
         plot.distr(data$`p.value`, main = paste('p-values -', metaname,type[f]), xlab = 'p-value')
      dev.off()

      # QQ plot.

      png(paste0(analyzedata[f], '_QC_Chi_lambda.png'))
         lambda <- qchisq(median(data$`p.value`), df = 1, lower.tail = FALSE) / qchisq(0.5, 1)
      dev.off()

      png(paste0(analyzedata[f], '_QC_QQplot_plot.png'))
         qqman::qq(data$`p.value`, main = sprintf('QQ plot of %s %s (lambda = %f)', metaname, type[f], lambda))
      dev.off()

      # Volcano plot.
      bt <- 0.01
      pt <- 3
      colors <- rep(gray(0.75, 0.25), nCpG)
      colors[data$beta >  bt & -log10(data$`p.value`) > 3] <- gg_colors(2, 0.5)[1] # Red
      colors[data$beta < -bt & -log10(data$`p.value`) > 3] <- gg_colors(2, 0.5)[2] # Blue

      png(paste0(analyzedata[f], '_QC_Volcano.png'))
         plot(data$beta, -log10(data$`p.value`), col = colors,
              main = paste('Volcano plot of', metaname), xlab = 'Beta', ylab = '-log10 p-value')
         abline(h = pt, v = c(-bt, bt), lty = 'dotted')
      dev.off()

      modeldata[[f]] <- data


   }

   names(modeldata) <- type;
   modeldata

}



# Annotate CpGs with Illumina 450K or EPIC data
annotate_CpGs <- function(CpGs, arrray)
{

   # Get library for 450K or Epic
   if( arrray == '450K'){
      library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
      ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
   }else{
      library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
      ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
   }

   # Download data

   ann$gene <- rownames(ann)

   # Merge cpgs with annotations
   CpGs.annot <-  merge( as.data.frame(CpGs), ann, by.x ="CpGs", by.y = "gene")

   # Remove data that we are not interested in (duplicate data column) before return
   ##..## return(CpGs.annot[ ,-(1)])
}
