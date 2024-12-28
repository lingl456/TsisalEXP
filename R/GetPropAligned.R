#' Title GetPropAligned
#'
#'#' Align Proportions Based on Correlation Coefficient
#'
#' This function aligns the columns of a given input matrix to a reference matrix based
#' on the correlation coefficients between their columns. The alignment is performed by
#' finding the maximum correlation for each column in the input matrix with columns in
#' the reference matrix. The function returns a transformed input matrix with columns
#' reordered according to their best correlation matches with the reference.
#'
#' This method is useful when you have a reference dataset, and you want to align or
#' reorder another dataset such that the columns of the input matrix match as closely
#' as possible to the reference matrix based on correlation.
#'
#' @param input A numeric matrix where the columns represent different samples or
#' features that need to be aligned with the reference matrix.
#' @param reference A numeric matrix that serves as the reference for alignment.
#' The columns of the reference matrix should correspond to those in the input matrix.
#'
#' @return A matrix with the columns of the input matrix reordered to align with
#' the reference matrix based on correlation.
#' @export
#'
#' @examples # Example of aligning an input matrix with a reference matrix based on correlation
#' input_data <- matrix(rnorm(100), nrow = 10, ncol = 5)  # Random data for illustration
#' reference_data <- matrix(rnorm(100), nrow = 10, ncol = 5)
#' aligned_data <- GetPropAligned(input_data, reference_data)
GetPropAligned <- function(input,reference){
  L = ncol(reference)
  colnames(input)=colnames(reference)=seq(1,dim(input)[2],by=1)
  corMat = cor(input,reference,use="pairwise.complete.obs")
  prop_cor = rep(0,L)
  tmpmat=corMat
  tmpmat[is.na(tmpmat)] = rnorm(L,-1,0.01)
  if(L>2){
    for(i in 1:L){
      maxind = which(tmpmat == max(tmpmat), arr.ind = TRUE)
      prop_cor[maxind[1]] = colnames(corMat)[maxind[2]]
      tmpmat[maxind[1],]=tmpmat[,maxind[2]]=rep(-1,L)
    }
  }else if(L==2){
    if(tmpmat[1,1]>0){
      prop_cor = c("1","2")
    }else{
      prop_cor = c("2","1")
    }
  }
  colnames(input) = prop_cor
  trans_input = input[,colnames(reference)]
  return(trans_input)
}
