#' Title
#'
#'#' Compute Akaike Information Criterion (AIC) for Deconvolution Model
#'
#' This function computes the Akaike Information Criterion (AIC) for a given deconvolution model
#' by evaluating the residual sum of squares (RSS) and the number of parameters in the model.
#' The AIC is used to compare models with different numbers of parameters and select the best-fitting model.
#' @param estProp A numeric matrix representing the estimated proportions of each cell type for each sample.
#' Rows correspond to genes, and columns correspond to samples.
#' @param Y.raw A numeric matrix representing the raw gene expression data, where rows correspond to genes and columns to samples.
#'
#' @return A numeric value representing the Akaike Information Criterion (AIC) for the model.
#' @export
#'
#' @examples # Example usage of compute_aic
#' estProp <- matrix(runif(100), nrow = 50, ncol = 2)  # Example estimated proportions (50 genes x 2 cell types)
#' Y_raw <- matrix(rnorm(1000), nrow = 50, ncol = 20)  # Example gene expression matrix (50 genes x 20 samples)
#' aic_value <- compute_aic(estProp, Y_raw)
#' print(aic_value)
compute_aic <- function(estProp,Y.raw){

  K = ncol(estProp)
  Nsample = dim(Y.raw)[2]
  idx = apply(estProp,2,function(x) sum(x) == 0)
  estProp[,idx] = matrix(runif(Nsample*sum(idx),0.0001,0.0002),Nsample,sum(idx))
  estProf <- t(mycsfit(estProp, t(Y.raw))$ghat)
  tmpmat <- estProf %*% t(estProp)
  rss = norm(Y.raw-tmpmat,type = "F")^2
  nSample = ncol(Y.raw) * nrow(Y.raw)
  nParam = K*(nrow(Y.raw)+ncol(Y.raw))
  aic = nSample*log(rss/nSample)+ 2*nParam + (2*nParam*(nParam+1))/(nSample-nParam-1)

  return(aic)
}
