% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/removeBatch2.R
\name{removeBatch2}
\alias{removeBatch2}
\title{Title removeBatch2}
\usage{
removeBatch2(dat.raw, pure.mean)
}
\arguments{
\item{dat.raw}{A numeric matrix of raw TCGA data, with rows representing genes and columns representing samples.}

\item{pure.mean}{A numeric matrix of reference (pure) data, with rows representing genes and columns representing samples.
The reference data should have the same gene names as the raw data.}
}
\value{
A list containing:
\describe{
\item{Y.raw}{A numeric matrix of the batch-corrected TCGA data.}
\item{pure}{A numeric matrix of the batch-corrected reference data.}
}
}
\description{
#' Remove Batch Effects from TCGA Data
}
\details{
This function removes batch effects from TCGA data by applying the ComBat method from the \code{sva} package.
It assumes that the raw data and reference data are matched by gene names.
}
\examples{
# Example usage:
result <- removeBatch2(dat.raw = tcga_data, pure.mean = reference_data)
batch_corrected_data <- result$Y.raw
corrected_reference_data <- result$pure
}
