#' Title CornerToEstProp
#'
#'#' Convert Corner Matrix to Estimated Proportions
#'
#' This function calculates the estimated proportions for each sample based on the corner matrix obtained
#' from a previous algorithm (e.g., SISAL). It uses non-negative least squares (NNLS) to estimate the proportions
#' of the components that best explain each sample.
#'
#' @param corner A numeric matrix of size N_sample x p, where N_sample is the number of samples and p is the
#' number of components. Each row represents a sample, and each column corresponds to a component (corner).
#'
#' @return A numeric matrix of size N_sample x p, where each row represents the estimated proportions for
#' each sample across the p components. The proportions are normalized to sum to 1 for each sample, and all
#' values are constrained between 0 and 1.
#' @export
#'
#' @examples # Example usage of CornerToEstProp
#' corner <- matrix(runif(100), nrow = 10, ncol = 10)  # Example corner matrix (10 samples x 10 components)
#' estProportions <- CornerToEstProp(corner)
CornerToEstProp <- function(corner){
  N_sample = dim(corner)[1]
  tmp <- nnls(corner, rep(1,N_sample))
  estProp <- diag(as.numeric(tmp$x),length(as.numeric(tmp$x))) %*% t(corner)
  estProp[estProp < 0] = 0
  estProp[estProp > 1] = 1
  return(t(estProp/colSums(estProp)))
}
