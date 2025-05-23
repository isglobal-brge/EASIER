% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getUniqueGenes.R
\name{getUniqueGenes}
\alias{getUniqueGenes}
\title{Get unique genes}
\usage{
getUniqueGenes(x, entrez = TRUE)
}
\arguments{
\item{x}{Annotated genes vector}

\item{entrez}{boolean. Optional if TRUE (default option), use org.Hs.egSYMBOL2EG to return entrez id}
}
\value{
Vector with unique genes
}
\description{
Get unique genes from annotated results
}
