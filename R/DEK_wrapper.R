#' Title DEK_wrapper
#'
#'#' Differential Expression with Tissue Proportions
#'
#' This function performs differential expression analysis between tissues based on their
#' proportions in the mixed tissue samples. It uses a wrapper around the `DEKTissue` function
#' to calculate tissue-specific differential expression, considering both single and multiple
#' tissue comparisons.
#' @param K The number of tissues in the dataset.
#' @param Y A numeric matrix or a `SummarizedExperiment` object representing the expression
#' matrix. Rows are features (genes), and columns are samples.
#' @param Prop A numeric matrix representing the tissue proportions for each sample. Rows are
#' samples, and columns are tissues.
#' @param design A numeric vector specifying the experimental design. If it contains zeros,
#' pairwise tissue comparisons will be made.
#' @param var.threshold A numeric value indicating the threshold for variance filtering (default is 0.1).
#'
#' @return A list of differential expression results. The list contains either:
#' - Tissue-specific differential expression results, if `design` is not all zeros.
#' - Pairwise tissue comparisons, if `design` is all zeros.
#' The results are returned as lists, with names corresponding to either tissue labels or
#' pairwise tissue comparisons.
#' @export
#'
#' @examples # Example of using DEK_wrapper to perform differential expression:
#' K <- 3
#' Y <- matrix(rnorm(1000), nrow = 100, ncol = 10)  # Example expression data
#' Prop <- matrix(runif(30), nrow = 10, ncol = 3)  # Example tissue proportions
#' design <- c(0, 1, 0)  # Example design vector (e.g., for pairwise comparisons)
#' DE_results <- DEK_wrapper(K, Y, Prop, design)
DEK_wrapper <- function(K, Y, Prop, design, var.threshold=0.1){
  DEK.res <- list()

  if(!all(design==0)){
    listname = rep("nn",K)
    for(k in 1:K){
      DEK.res[[k]] = DEKTissue(K, Y, Prop, design, var.threshold=var.threshold, WhichPar=k+K)
      if(is.null(colnames(Prop))){
        listname[k] = paste0("Tissue",k)
      }else if(!is.null(colnames(Prop))){
        listname[k] = colnames(Prop)[k]
      }
    }

  }

  else if(all(design==0)){
    g=1
    listname = rep("nn",K*(K-1)/2)
    for(i in 1:(K-1)){
      for(j in 2:K){
        cc = rep(0,K)
        cc[i] = 1
        cc[j] = -1
        DEK.res[[g]] = DEKTissue(K, Y, Prop, design, var.threshold = var.threshold, contrast.vec = cc)
        listname[g] = paste0("Tissue",i,"vs",j)
        g = g+1
      }
    }
  }
  names(DEK.res) = listname
  return(DEK.res)
}
