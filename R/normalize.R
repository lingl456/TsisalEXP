#' Title normalize
#'
#'#' Normalize Data Matrix by Row Sums
#'
#' This function normalizes the input data matrix by dividing each element by the sum of its corresponding row.
#' This is useful for standardizing gene expression data, where the rows represent genes and columns represent samples.
#'
#' @param Y A numeric matrix of size L x N, where L is the number of genes and N is the number of samples.
#' Each element represents the expression value of a gene in a sample.
#'
#' @return A numeric matrix of the same size as Y, where each element is normalized by dividing it by the sum of its row.
#' @export
#'
#' @examples # Example usage of the normalize function
#' data <- matrix(rnorm(1000), nrow = 100, ncol = 10) # Random example data (100 genes x 10 samples)
#' normalized_data <- normalize(data)
normalize <- function(Y){
  Y.norm = Y/rowSums(Y)
  return(Y.norm)
}
