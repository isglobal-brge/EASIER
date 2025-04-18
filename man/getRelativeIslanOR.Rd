% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getRelativeIslandOR.R
\name{getRelativeIslanOR}
\alias{getRelativeIslanOR}
\title{Significative variable OR by position Relative to Island}
\usage{
getRelativeIslanOR(significative, position, varname)
}
\arguments{
\item{significative}{numerical. Variable to take in to account as significative variable, could be FDR, p-value,...}

\item{position}{vector with position Relative to Island}

\item{varname}{string. position Relative to Island name. For example : "Island", "OpenSea", "S_Shore", "N_Shore", "S_Shelf", "N_Shelf"}
}
\description{
Get OR by Position Relative to Island taking in to account a significative variable parameter like FDR values
}
