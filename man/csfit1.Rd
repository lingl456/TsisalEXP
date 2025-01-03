% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/csfit1.R
\name{csfit1}
\alias{csfit1}
\title{Title csfit1}
\usage{
csfit1(cc, G, logRm = FALSE, logBase = 2)
}
\arguments{
\item{cc}{A numeric matrix or vector representing the data (e.g., response variable). Each column corresponds
to a sample and each row corresponds to a feature.}

\item{G}{A numeric matrix or vector representing the predictor variable(s). It should have the same number of rows
as \code{cc}. Each column represents a different predictor.}

\item{logRm}{A logical value (default is FALSE). If TRUE, the function will apply a logarithmic transformation
to the matrix \code{G} before fitting the linear model.}

\item{logBase}{The base of the logarithm used for the transformation (default is 2). Only used if \code{logRm} is TRUE.}
}
\value{
A list containing two components:
\item{ghat}{A numeric vector or matrix of fitted coefficients (beta values) for the linear model. If \code{logRm} is TRUE,
the coefficients are back-transformed using the specified log base.}
\item{residuals}{A numeric vector or matrix of residuals from the linear regression fit.}
}
\description{
#' Perform Linear Regression for a Given Data Matrix
}
\details{
This function fits a linear model to the data, with an optional log transformation applied to the input
matrix before fitting. The model does not include an intercept term.
}
\examples{
# Generate random data
set.seed(123)
cc <- matrix(rnorm(100), nrow = 10, ncol = 10)
G <- matrix(rnorm(100), nrow = 10, ncol = 10)

# Fit linear model without log transformation
result <- csfit1(cc, G)

# Fit linear model with log transformation (base 2)
result_log <- csfit1(cc, G, logRm = TRUE, logBase = 2)
}
