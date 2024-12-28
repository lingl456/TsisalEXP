#' Title mysisal
#'
#' SISAL-based Tissue Proportion Estimation and Marker Selection
#'
#' This function uses SISAL (Sparse Inverse Simultaneous Component Analysis) to estimate tissue
#' proportions from gene expression data and to identify the top N markers for each tissue type.
#'
#' @param Y  A G*N matrix where G is the number of genes, and N is the number of samples. This matrix
#' represents gene expression data.
#' @param K  A G*N matrix where G is the number of genes, and N is the number of samples. This matrix
#' represents gene expression data.
#' @param topN The number of top markers to select based on their contribution to the tissue proportions.
#'
#' @return A list containing:
#' - `estProp`: A matrix of estimated tissue proportions for each sample.
#' - `selMarker`: A vector of indices corresponding to the top N markers selected based on the SISAL method.
#' @export
#'
#' @examples # Example usage of mysisal function:
#' Y <- matrix(rnorm(1000), nrow = 100, ncol = 10)  # Example gene expression data
#' K <- 3  # Assume there are 3 tissue types
#' topN <- 10  # Select top 10 markers
#' result <- mysisal(Y, K, topN)
#'
#' # Access the estimated tissue proportions and selected markers
#' result$estProp
#' result$selMarker
mysisal <- function(Y, K, topN){
  Y.norm = normalize(Y)
  sisalres = sisal(t(Y.norm), p = K, iters = 100)
  corner = sisalres$endpoints
  distances = sisalres$distances
  estProp = CornerToEstProp(corner)
  selMarker = CornerToMarker(distances, topN)
  return(list(estProp = estProp,selMarker = selMarker))
}
