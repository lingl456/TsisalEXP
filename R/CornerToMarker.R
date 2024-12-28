#' Title
#'
#'#' Identify Marker Genes Based on Distances
#'
#' This function identifies the top N marker genes for each component based on their distances in a given distance matrix.
#' For each component, the function selects the top N genes with the smallest distance (i.e., most closely related) and
#' considers them as marker genes.
#'
#' @param distances A numeric matrix of size N_genes x p, where N_genes is the number of genes and p is the number of components.
#' Each element represents the distance between a gene and a component. Smaller values indicate higher similarity between the gene and the component.
#' @param topN An integer specifying the number of marker genes to select for each component.
#'
#' @return A list of length p (number of components), where each element is a character vector containing the names of the
#' top N marker genes for the corresponding component. The marker genes are selected based on the smallest distances for each component.
#' @export
#'
#' @examples # Example usage of CornerToMarker
#' distances <- matrix(runif(100), nrow = 10, ncol = 10)  # Example distance matrix (10 genes x 10 components)
#' topN <- 3
#' markerGenes <- CornerToMarker(distances, topN)
CornerToMarker <- function(distances,topN){
  markerList <- apply(distances, 2, function(xx) {
    pure <- rownames(distances)[order(xx)[1:topN]]
    return(pure)
  })
  markerList <-  split(markerList, rep(1:ncol(markerList), each = nrow(markerList)))
  return(markerList)
}
