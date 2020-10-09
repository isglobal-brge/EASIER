#' Annotate a List of CpGs with Illumina 450K or EPIC data
#'
#' Get CpG annotations from Illumina 450K or EPIC array type
#'
#' @param CpGs string. Path with results from GWAMA
#' @param artype string, Illumina array type, 450K or EPIC
#' @param filename string vector with file name where annotations
#' @param outdir string vector. Route to output folder to store results
#'
#' @return dataframe with annotated CpGs
#'
#' @export
get_annotattions <- function(CpGs, artype, filename, outdir)
{

   # Get library for 450K or Epic
   if( artype == '450K')
      library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
   else
      library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

   # Output filename
   outfilename <- tools::file_path_sans_ext(basename(filename))
   outfilename <- paste(outdir, outfilename, sep="/")

   # Download data
   ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
   ann$gene <- rownames(ann)

   # Merge cpgs with annotations
   CpGs.annot <-  merge( as.data.frame(CpGs), ann, by.x ="CpGs", by.y = "gene")

   # Remove data that we are not interested in (duplicate data) before write results
   ##..## CpGs.annot <- CpGs.annot[ ,-(1)]

   # Write data to a file
   write.table( CpGs.annot, paste0(outfilename,"_annotated.txt"), quote=F, sep="\t")

   # Return data
   return(CpGs.annot)
}
