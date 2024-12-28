#' Title getProportion
#'
#' Generate Proportions for Dirichlet Distribution
#'
#' This function generates proportions for a given number of samples based on a Dirichlet distribution.
#' The Dirichlet distribution is parameterized using pre-defined concentration parameters (`alpha.ctr`) that are
#' either fixed for 4 components or can be adjusted for a different number of components (`L`).
#' The generated proportions are typically used to model the relative abundance of different categories or tissues.
#'
#' @param N_sample An integer specifying the number of samples to generate. Each sample will have a proportion vector.
#' @param L An integer specifying the number of categories or components. The default is 4.
#'
#' @return A numeric matrix of size `N_sample x L`, where each row represents a sample and the columns represent the proportions
#' for each category/component. The values in each row sum to 1, as they represent proportions.
#'
#' @export
#'
#' @examples #' # Generate proportions for 100 samples with 4 categories
#' proportions <- getProportion(100)
#' head(proportions)
#'
#' # Generate proportions for 50 samples with 6 categories
#' proportions_6 <- getProportion(50, L = 6)
#' head(proportions_6)
getProportion <- function(N_sample, L=4){

  # get MLE for dirichlet distribution
  if(L == 4) {
    alpha.ctr=c(0.968, 4.706, 0.496, 0.347)
  } else {
    alpha.ctr=c(0.89, 4.12, 0.47, 0.33, 0.61, 1.02, 2.25)[1:L]
  }

  prop.matrix.ctr = rdirichlet(N_sample, alpha.ctr)

  return(prop.matrix.ctr)
}
