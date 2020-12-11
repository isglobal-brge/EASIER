#' Venn diagrams plots
#'
#' Plot Venn diagrams
#'
#' @param venndata prefixes used as a file name to store QC data
#' @param path optional, path where QC data is stored, by default current path
#' @param qcpath optional, path to store Venn Diagrams, by default current path
#' @param pattern optional, pattern with data to be used in Venn diagrams, by default '*_QCData.txt'
#' @param bn, optional, column name with bonferroni adjusted values, if NULL, no venndiagram plot for Bonferroni
#' @param fdr, column name with bonferroni adjusted values, if NULL, no venndiagram plot for FDR
#' @param threshold, numeric, optiona, by default threshlold = 0.05, we take in to account all CpGs with p-value lower than this threshold
#'
#' @return Venn diagrams
#'
#' @import VennDiagram RColorBrewer
#'
#' @export
plot_venndiagram <- function(venndata, qcpath = '.', plotpath = '.', pattern = '*_QCData.txt',  bn = NULL, fdr = NULL, threshold = 0.05)
{
   if( is.null(bn)  & is.null(fdr) )
      stop("No Venn diagram to plot")

   if(is.null(venndata))
      stop("No data to plot")

   listnames <- venndata
   myCol <- RColorBrewer::brewer.pal(length(listnames), "Set3")

   cohortsdata <- multmerge(qcpath, venndata, pattern)

   if(! bn %in% colnames(cohortsdata[[1]]))
      warning("Field for Bonferroni not found in data")

   if(!fdr %in% colnames(cohortsdata[[1]]))
      warning("Field for FDR not found in data")

   if(!is.null(bn) )
   {
      selection <- paste0( "lapply(cohortsdata, subset,", paste0('`',bn,'`')," =='yes')")

      # Venn Diagram - Bonferroni
      VennDiagram::venn.diagram(
         x = lapply(lapply( eval(parse(text = selection))  , "[", 1),unlist),
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

   if( !is.null(fdr) )
   {
      selection <- paste0( "lapply(cohortsdata, subset,", paste0('`',fdr,'`')," <",threshold,")")

      VennDiagram::venn.diagram(
         x = lapply(lapply(eval(parse(text = selection)), "[", 1),unlist),
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
