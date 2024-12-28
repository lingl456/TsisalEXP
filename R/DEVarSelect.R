#' Title DEVarSelect
#'
#'#' Select Differentially Expressed Marker Genes for Deconvolution
#'
#' This function selects marker genes based on their differential expression across cell types.
#' It uses a statistical method to identify tissue-specific genes, which are then used in
#' the deconvolution process. The function selects a set of markers for each cell type
#' to improve the estimation of cell type proportions in mixed tissue samples.
#'
#' @param Y_raw A numeric matrix or a `SummarizedExperiment` object. The rows represent genes
#' and the columns represent samples. This is the input data from complex tissues for marker
#' gene selection.
#' @param Prop0 A numeric matrix of proportions, with rows representing samples and columns
#' representing cell types. These are the initial estimates of the cell type proportions.
#' @param nMarker The number of marker genes to select for each cell type. Default is 1000.
#' @param bound_negative A logical value indicating whether to bound the selected proportions
#' to be non-negative. Default is `FALSE`.
#'
#' @return A vector of indices corresponding to the selected marker genes. These genes are
#' expected to be the most differential across the tissue types and useful for the deconvolution process.
#' @export
#'
#' @examples # Example of selecting marker genes for deconvolution:
#' Y_raw <- matrix(rnorm(1000), nrow = 100, ncol = 10)  # Example data
#' Prop0 <- matrix(runif(30), nrow = 10, ncol = 3)  # Example cell type proportions
#' marker_genes <- DEVarSelect(Y_raw, Prop0, nMarker = 500)
#'
#' # Access the selected marker genes
#' marker_genes
DEVarSelect <- function(Y_raw, Prop0, nMarker = 1000, bound_negative = FALSE){

  if (is(Y_raw, "SummarizedExperiment")) {
    se <- Y_raw
    Y_raw <- assays(se)$counts
  } else if (!is(Y_raw, "matrix")) {
    stop("Y_raw should be a matrix or a SummarizedExperiment object!")
  }

  if (nrow(Prop0) < ncol(Prop0)) {
    stop("Prop0 should have dimension N (samples) by K (cell types)!")
  }
  if (!ncol(Y_raw) == nrow(Prop0)) {
    stop("Y_raw should have dimension P (features) by N (samples)!")
  }

  K <- ncol(Prop0)
  N_sample <- nrow(Prop0)

  ## find tissue specific genes
  idx <- NULL
  for(k in seq_len(K)) {
    cvec <- rep(-1/(K-1),K)
    cvec[k] <- 1
    design <- rep(0,N_sample)
    tmp <- DEKTissue(K, Y=Y_raw,
                     Prop=Prop0,
                     design=design,
                     contrast.vec=cvec)
    idx[[k]] <- sort(abs(tmp$t.stat),
                     decreasing=TRUE,
                     index=TRUE)$ix
  }
  nmarker <- nMarker
  ## number of markers per tissue. Consider overlaps
  nmarker_tissue <- nmarker/K * 1.2
  idxMarker <- NULL
  for(k in seq_len(K)) {
    idxMarker <- c(idxMarker,
                   idx[[k]][seq_len(nmarker_tissue)])
  }
  idxMarker <- unique(idxMarker)

  return(idxMarker)
}
