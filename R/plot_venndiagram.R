#' Venn diagrams plots
#'
#' Plot Venn diagrams
#'
#' @param venndata prefixes used as a file name to store QC data
#' @param path optional, path where QC data is stored, by default current path
#' @param qcpath optional, path to store Venn Diagrams, by default current path
#' @param bn, optional, if true gets venn diagram with bonferroni adjusted values
#' @param fdr, optional, if true gets venn diagram with fdr adjusted values
#'
#' @return Venn diagrams
#'
#' @export
plot_venndiagram <- function(venndata, qcpath = '.', plotpath = '.', bn = TRUE, fdr = TRUE)
{
   if(bn==FALSE & fdr == FALSE)
      stop("No Venn diagram to plot")

   if(is.null(venndata))
      stop("No data to plot")

   listnames <- venndata
   myCol <- RColorBrewer::brewer.pal(length(listnames), "Set3")

   # Get all data and put it to a list
   filenames <- list.files(path = qcpath, pattern = '*_QCData.txt', full.names = TRUE)

   if(is.null(filenames))
      stop( paste("No QC data found in",qcpath))

   cohortsdata <- multmerge(qcpath, venndata)

   if(bn == TRUE)
   {
      # Venn Diagram - Bonferroni
      VennDiagram::venn.diagram(
         x = lapply(lapply(lapply(cohortsdata, subset, `padj.bonf`=='yes'), "[", 1),unlist),

         # numbers
         fontfamily = "symbol",

         # Cattegory
         cat.fontfamily = "sans", cat.cex = 0.3,
         cat.default.pos = "outer",
         category.names = listnames,

         # Circles
         lwd = 2,  fill = myCol,

         # Output features
         height = 1200, width = 1200, resolution = 300, imagetype = 'png',

         filename = paste0(plotpath,'/venn_diagramm_bonfer_',i,'.png'),
         output=TRUE
      )
   }

   # Venn Diagram - FDR

   if(fdr == TRUE)
   {
      VennDiagram::venn.diagram(
         x = lapply(lapply(lapply(cohortsdata, subset, `padj.fdr`==0.05), "[", 1),unlist),

         # numbers
         fontfamily = "symbol",

         # Cattegory
         cat.fontfamily = "sans", cat.cex = 0.3,
         cat.default.pos = "outer",
         category.names = listnames,

         # Circles
         lwd = 2,  fill = myCol,

         # Output features
         height = 1200, width = 1200, resolution = 300, imagetype = 'png',

         filename = paste0(plotpath,'/venn_diagramm_fdr_',i,'.png'),
         output=TRUE
      )
   }



}


# Function to read all files with certain "pattern" in specified "path" and put it to a list
multmerge <-  function(mypath, arfiles)
{
   pattern <- '_QCData.txt'
   filenames <- paste(mypath,paste0(arfiles,pattern), sep="/")
   datalist = lapply(filenames, function(x){read.delim2(file=x, header=TRUE, stringsAsFactors = FALSE)})
   return(datalist)
}
