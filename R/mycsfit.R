#' Title
#'
#'#' Fit a Linear Model to Cell Type Proportions
#'
#' This function fits a linear model to cell type proportions (`cc`) and gene expression data (`G`),
#' using least squares regression. The function allows for optional log-transformation of the gene expression data
#' and provides the fitted coefficients, residuals, standard errors, and other related statistics.
#'
#' @param cc A numeric matrix or data frame representing cell type proportions for each sample.
#' Rows correspond to cell types, and columns correspond to samples.
#' @param G A numeric matrix or data frame representing the gene expression data, where rows correspond to genes and columns to samples.
#' @param logRm A logical value indicating whether to apply a log transformation to the gene expression data. Default is `FALSE`.
#' @param logBase The base of the logarithm to use if `logRm` is `TRUE`. Default is 2 (log base 2).
#'
#' @return A list of class `csfit` containing the following components:
#'   - `ghat`: A numeric matrix of the fitted coefficients for each cell type and gene.
#'   - `residuals`: A numeric matrix of the residuals from the regression for each gene and sample.
#'   - `se`: A numeric matrix of the standard errors for the fitted coefficients.
#'   - `n`: The number of samples in the gene expression data `G`.
#' @export
#'
#' @examples # Example usage of mycsfit
#' cc <- matrix(runif(100), nrow = 10, ncol = 10)  # Example cell type proportions (10 cell types x 10 samples)
#' G <- matrix(rnorm(1000), nrow = 100, ncol = 10)  # Example gene expression matrix (100 genes x 10 samples)
#' result <- mycsfit(cc, G)
#' print(result$ghat)
mycsfit <- function(cc,G,logRm=FALSE,logBase=2) {
  if(logRm == TRUE) {
    G <- logBase^G
  }

  fit1 <- lsfit(cc, G, intercept = FALSE)
  se1 <- ls.diag(fit1)$std.err
  if(logRm == TRUE) {
    ghat <- log(coef(fit1),logBase)
    ghat[is.nan(ghat)] <- 0
    se <- log(se1,logBase)
    res <- list(ghat = ghat, residuals = residuals(fit1),se = se1, n = nrow(G))
  }
  else {
    res <- list(ghat = coef(fit1), residuals = residuals(fit1),se= se1, n = nrow(G))
  }

  # fix dimnames
  dimnames(res$residuals) <- dimnames(G)
  features <- colnames(G)
  colnames(res$ghat) <- features
  colnames(res$se) <- features

  # return
  class(res) <- 'csfit'
  res
}

