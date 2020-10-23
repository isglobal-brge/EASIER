#' Run GWAMA meta-analysis
#'
#' Runs GWAMA for fixed and random effects, results with fixed events are stored in file wit _Fixed suffix and meta-analysis results with random effects are in a sufixed file _Random.
#'
#' @param inputfolder string. Path with input files for GWAMA.
#' @param outputfolder string. Path to store gwama meta-analisis results
#' @param outputfilename string. File name to write gwama meta-analisis results, will end with _Fixed or with _Random.
#' @param gwama.dir string. Route to GWAMA binary
#' @param hapmapfile string. complete or relative route to hapmap file
#'
#' @return Route to gwama result files
#'
#'
#' @export
run_GWAMA_MetaAnalysis <- function(inputfolder, outputfolder, outputfilename, gwama.dir)
{

   # -- GWAMA meta-analysis --

   outputfolder <- ifelse( substr(outputfolder, 1, 2) == './', substr(outputfolder,3,nchar(outputfolder)), outputfolder)

   # Check if gwama results folder exists and create it.
   if(!dir.exists(file.path(getwd(), outputfolder ))) suppressWarnings(dir.create(file.path(getwd(), outputfolder), recursive = TRUE))

   # Fix effects
   # gwm.fixedoufile <- paste0(str_replace(getwd(),' ','\\\\ '), "/", outputfolder, "/", outputfilename,"_Fixed" )
   gwm.fixedoufile <- paste0(getwd(), "/", outputfolder, "/", outputfilename,"_Fixed" )
   # gwamaexec.fix <- paste0( gwama.dir, "GWAMA -qt -m ", hapmapfile,
   #                          " -i ", paste(str_replace(getwd(),' ','\\\\ '), prefixgwama, "gwama.in", sep="/"),
   #                          " -o ", gwm.fixedoufile)

   gwamaexec.fix <- paste0( gwama.dir, "GWAMA -qt -m ", hapmapfile,
                            " -i ", paste(getwd(), inputfolder, "gwama.in", sep="/"),
                            " -o ", gwm.fixedoufile)

   gwamaexec.fix
   system(gwamaexec.fix)

   # Random effects
   gwm.ramdomoufile <- paste0(getwd(), "/", outputfolder, "/", outputfilename,"_Random" )
   # gwamaexec.random <- paste0( gwama.dir, "GWAMA -r -qt -m ", hapmapfile,
   #                             " -i ", paste(getwd(),prefixgwama, "gwama.in", sep="/"),
   #                             " -o ", gwm.ramdomoufile)

   gwamaexec.random <- paste0( gwama.dir, "GWAMA -r -qt -m ", hapmapfile,
                               " -i ", paste(getwd(), inputfolder, "gwama.in", sep="/"),
                               " -o ", gwm.ramdomoufile)
   system(gwamaexec.random)


   # -- GWAMA PLOTS : Manhattan and QQ-Plots  --

   qqplot <- paste0(" < ", gwama.dir, "QQ.R")
   manhattanplot <- paste0(" < ", gwama.dir, "MANH.R")

   # Exec QQ-plot for Fixed and Random effects
   plotexec.fix <- paste0( "R --slave --vanilla --args input=", paste(getwd(), outputfolder, paste0(outputfilename,"_Fixed.out"),sep="/" ), " out=", paste0(getwd(), "/", outputfolder, "/", outputfilename,"_Fixed.QQ.png" ), qqplot)
   system(plotexec.fix)
   plotexec.random <- paste0( "R --slave --vanilla --args input=", paste(getwd(), outputfolder, paste0(outputfilename,"_Random.out"),sep="/" ), " out=", paste0(getwd(), "/", outputfolder, "/", outputfilename,"_Random.QQ.png" ), qqplot)
   system(plotexec.random)

   # Exec Manhattan-plot for Fixed and Random effects
   plotexec.fix <- paste0( "R --slave --vanilla --args input=", paste(getwd(), outputfolder, paste0(outputfilename,"_Fixed.out"),sep="/" ), " out=", paste0(getwd(), "/", outputfolder, "/", outputfilename,"_Fixed.Manhattah.png" ), manhattanplot)
   system(plotexec.fix)
   plotexec.random <- paste0( "R --slave --vanilla --args input=", paste(getwd(), outputfolder, paste0(outputfilename,"_Random.out"),sep="/" ), " out=", paste0(getwd(), "/", outputfolder, "/", outputfilename,"_Random.Manhattah.png" ), manhattanplot)
   system(plotexec.random)

   files <- c(gwm.fixedoufile, gwm.ramdomoufile)
   return(files)
}

