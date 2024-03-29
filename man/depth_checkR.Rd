% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qualityflag.R
\name{depth_checkR}
\alias{depth_checkR}
\title{Checks stability of sample by comparison of depth to vpos measurements on sonde.}
\usage{
depth_checkR(x)
}
\arguments{
\item{x}{A tibble as exported from either \code{\link{imp_raw_csv}}
or \code{\link{old_data_updatR}}}
}
\value{
A tibble where \code{QUAL} variable has 1000 added if > 0.1m variation.
}
\description{
Checks stability of sample by comparison of depth to vpos measurements on sonde.
}
