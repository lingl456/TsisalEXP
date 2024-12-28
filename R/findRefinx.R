#' Title findRefinx
#'
#'#' Select Marker Features Based on Variance or Coefficient of Variation
#'
#' This function selects marker features (genes) from a given dataset based on either their
#' coefficient of variation (CV) or variance. The markers are selected based on their
#' variability across samples, which can be useful for downstream analysis such as
#' dimensionality reduction or feature selection.
#'
#' The function supports both raw count data from RNA sequencing or other types of genomic data
#' stored as either a matrix or a `SummarizedExperiment` object. The features are sorted based on
#' their variability, and the top `nmarker` features are returned.
#'
#' @param rawdata A matrix or `SummarizedExperiment` object containing the raw data.
#' The rows represent features (genes), and the columns represent samples.
#' If a `SummarizedExperiment` is provided, the raw count data is retrieved from the assay.
#' @param nmarker An integer specifying the number of marker features to select based on their variability.
#' Defaults to 1000.
#' @param sortBy A character string specifying the criterion to sort the features by.
#' Options are "cv" for coefficient of variation or "var" for variance. Defaults to "var".
#'
#'
#' @return A vector of indices corresponding to the top `nmarker` most variable features,
#' sorted according to the selected criterion.
#' @export
#'
#' @examples # Example of selecting top 1000 most variable features from a count matrix
#' rawdata <- matrix(rnorm(10000), nrow = 100, ncol = 100)  # Random data for illustration
#' top_features <- findRefinx(rawdata, nmarker = 1000, sortBy = "var")
#' top_features  # View the indices of selected features
findRefinx <- function(rawdata,
                       nmarker = 1000,
                       sortBy = "var") {

  if (is(rawdata, "SummarizedExperiment")) {
    se <- rawdata
    rawdata <- assays(se)$counts
  } else if (!is(rawdata, "matrix")) {
    stop("rawdata should be a matrix
               or a SummarizedExperiment object!")
  }

  if (nrow(rawdata) < ncol(rawdata)) {
    stop("rawdata matrix should have
               dimension P (features) by N (samples)!")
  }
  if (nmarker > nrow(rawdata)) {
    stop("You have specified nmarker larger
               than the number of original markers!")
  }

  if (sortBy == "cv") {
    mm <- Matrix::rowMeans(rawdata)
    vv <- matrixStats::rowVars(rawdata)
    cv <- sqrt(vv) / mm
    cv[is.na(cv)] <- 0
    final_v <- cv
  } else if(sortBy == "var") {
    vv <- matrixStats::rowVars(log(rawdata+1))
    vv[is.na(vv)] <- 0
    final_v <- vv
  } else {
    stop("sortBy should be either 'cv' or 'var'!")
  }

  ix <- sort(final_v, decreasing=TRUE, index=TRUE)$ix
  return(ix[seq(nmarker)])
}
