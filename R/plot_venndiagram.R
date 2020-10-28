#' Venn diagrams plots
#'
#' Plot Venn diagrams
#'
#' @param venndata prefixes used as a file name to store QC data
#' @param path optional, path where QC data is stored, by default current path
#' @param qcpath optional, path to store Venn Diagrams, by default current path
#' @param pattern optional, pattern with data to be used in Venn diagrams, by default '*_QCData.txt'
#' @param bn, optional, if true gets venn diagram with bonferroni adjusted values
#' @param fdr, optional, if true gets venn diagram with fdr adjusted values
#'
#' @return Venn diagrams
#'
#' @export
plot_venndiagram <- function(venndata, qcpath = '.', plotpath = '.', pattern = '*_QCData.txt', bn = TRUE, fdr = TRUE)
{
   if(bn==FALSE & fdr == FALSE)
      stop("No Venn diagram to plot")

   if(is.null(venndata))
      stop("No data to plot")

   listnames <- venndata
   myCol <- RColorBrewer::brewer.pal(length(listnames), "Set3")

   cohortsdata <- multmerge(qcpath, venndata, pattern)

   if(bn == TRUE)
   {
      # Venn Diagram - Bonferroni
      VennDiagram::venn.diagram(
         x = lapply(lapply(lapply(cohortsdata, subset, `padj.bonf`=='yes'), "[", 1),unlist),
         filename = paste0(plotpath,'/venn_diagramm_bonfer_',i,'.png'),
         output=TRUE,

         # Circles
         lwd = 2,  fill = myCol,

         # # numbers
         cex = .6,
         fontfamily = "sans",

         # Cattegory
         cat.fontfamily = "sans", cat.cex = 0.3,
         cat.default.pos = "outer",
         category.names = listnames,

         # Output features
         height = 1200, width = 1200, resolution = 300, imagetype = 'png'
      )
   }

   # Venn Diagram - FDR

   if(fdr == TRUE)
   {
      VennDiagram::venn.diagram(
         x = lapply(lapply(lapply(cohortsdata, subset, `padj.fdr`<0.05), "[", 1),unlist),
         filename = paste0(plotpath,'/venn_diagramm_fdr_',i,'.png'),
         output=TRUE,

         # Circles
         lwd = 2,  fill = myCol,

         # numbers
         cex = .6,
         fontfamily = "sans",

         # Cattegory
         cat.fontfamily = "sans", cat.cex = 0.3,
         cat.default.pos = "outer",
         category.names = listnames,

         # Output features
         height = 1200, width = 1200, resolution = 300, imagetype = 'png'
      )
   }



}


# Function to read all files with certain "pattern" in specified "path" and put it to a list
multmerge <-  function(mypath, arfiles, patt)
{

   filenames <- list()
   filenames <- list.files(path = mypath, pattern = patt, full.names = TRUE)

   if(length(filenames)==0)
      filenames <- unlist(lapply(arfiles, function(a) {list.files(file.path(mypath,a), pattern = patt, full.names = TRUE) }))

   if(!exists("filenames") & length(filenames)>0)
      stop( paste("No QC data found in",mypath))

   datalist = lapply(filenames, function(x){read.delim2(file=x, header=TRUE, stringsAsFactors = FALSE)})

   return(datalist)
}
