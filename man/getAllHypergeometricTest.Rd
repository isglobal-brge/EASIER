% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getAllHypergeometricTest.R
\name{getAllHypergeometricTest}
\alias{getAllHypergeometricTest}
\title{HyperGeometric Test for positions relative to Island}
\usage{
getAllHypergeometricTest(
  significative,
  criteria,
  outputdir = ".",
  outputfile = NULL,
  plots = TRUE
)
}
\arguments{
\item{significative}{numerical. Variable to take in to account as significative variable, could be FDR, p-value,...}

\item{criteria}{dataframe with all relative positions to island values or chromatin status 15 or 18 Fetal Placenta}

\item{outputdir}{string. Output path to store file with results, by default results are written in current dir}

\item{outputfile}{string. File name to store results if no name is provided results are not written.}

\item{plots}{boolean. If plot is TRUE, plot results}
}
\value{
Dataframe with hypergeometric test results for all criteria values
}
\description{
HyperGeometric Test for positions relative to Island
}
