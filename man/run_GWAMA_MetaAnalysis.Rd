% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_GWAMA_MetaAnalysis.R
\name{run_GWAMA_MetaAnalysis}
\alias{run_GWAMA_MetaAnalysis}
\title{Run GWAMA meta-analysis}
\usage{
run_GWAMA_MetaAnalysis(
  inputfolder,
  outputfolder,
  outputfilename,
  gwama.dir,
  hapmapfile
)
}
\arguments{
\item{inputfolder}{string. Path with input files for GWAMA.}

\item{outputfolder}{string. Path to store gwama meta-analisis results}

\item{outputfilename}{string. File name to write gwama meta-analisis results, will end with _Fixed or with _Random.}

\item{gwama.dir}{string. Route to GWAMA binary}

\item{hapmapfile}{string. complete or relative route to hapmap file}
}
\value{
Route to gwama result files
}
\description{
Runs GWAMA for fixed and random effects, results with fixed events are stored in file wit _Fixed suffix and meta-analysis results with random effects are in a sufixed file _Random.
}
