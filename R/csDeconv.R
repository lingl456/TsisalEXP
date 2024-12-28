#' Title csDeconv
#'
#'#' Cell Type Deconvolution Using Reference-Free Deconvolution Method
#'
#' This function performs cell type deconvolution on high-throughput tissue measurement data.
#' The deconvolution process aims to estimate the proportion of pure cell types present in
#' each sample based on gene expression data from complex tissue samples. The deconvolution
#' is reference-free, and the function supports iterative refinement to improve proportion
#' estimates.
#'
#' @param Y_raw A numeric matrix or a `SummarizedExperiment` object. The rows represent genes
#' and the columns represent samples. This is the input data from complex tissues for  deconvolution.
#' @param K The number of pure cell types to be estimated from the data.
#' @param nMarker The number of marker genes to use for the deconvolution. Default is 1000.
#' @param InitMarker A vector of marker gene names to be used for the deconvolution. If not
#' provided, the top variable genes will be selected based on the input data.
#' @param TotalIter The total number of iterations to perform for the deconvolution process.
#' Default is 30.
#' @param bound_negative The total number of iterations to perform for the deconvolution process.
#' Default is 30.
#'
#' @return A list containing the following components:
#' \item{allRMSE}{A vector of RMSE values for each iteration of the deconvolution process.}
#' \item{allProp}{A list containing the estimated proportions for each iteration.}
#' \item{estProp}{The estimated proportions of pure cell types for the sample corresponding
#' to the iteration with the minimum RMSE.}
#' \item{updatedInx}{The indices of the marker genes used for the deconvolution.}
#'
#' @export
#'
#' @examples # Example of deconvolution with random data:
#' Y_raw <- matrix(rnorm(1000), nrow = 100, ncol = 10)  # Example data
#' K <- 3  # Number of cell types
#' deconv_result <- csDeconv(Y_raw, K, nMarker = 500)
#'
#' # Access the estimated proportions
#' deconv_result$estProp
csDeconv <- function(Y_raw,
                     K,
                     nMarker = 1000,
                     InitMarker = NULL,
                     TotalIter = 30,
                     bound_negative = FALSE) {
  # Y_raw is the high-throughput measurement from complex tissues
  #      (rows for features and columns for samples);
  #      or a SummarizedExperiment object.
  # K is the pre-specified number of pure cell types
  # FUN is the reference-free deconvolution function,
  #    this function should take Y_raw and K,
  #    and the return values should be a N by K proportion matrix.
  #    N is the number of samples and K is the number of cell types.
  # nMarker is the number of marker used in the deconvolution
  # InitMarker is the initial marker used in the deconvolution,
  #          if not specified, the top variable features will be used
  # TotalIter is the total number of iterations specified

  if (is(Y_raw, "SummarizedExperiment")) {
    se <- Y_raw
    Y_raw <- assays(se)$counts
  } else if (!is(Y_raw, "matrix")) {
    stop("Y_raw should be a matrix or a SummarizedExperiment object!")
  }

  if (is.null(rownames(Y_raw))) {
    row.names(Y_raw) <- seq(nrow(Y_raw))
  }
  if (is.null(InitMarker)) {
    if (nrow(Y_raw) < 2000) {
      InitMarker <- findRefinx(Y_raw, nmarker = nMarker)
    } else {
      tmp <- findRefinx(Y_raw, nmarker = nMarker*2)
      InitMarker <- tmp[nMarker+1:nMarker]
    }
  } else {
    if (sum(!(InitMarker %in% rownames(Y_raw))) > 0) {
      stop("Discrepancy between
           InitMarker and the row names of Y_raw!")
    }
  }

  allProp <- list()
  allRMSE <- rep(0, TotalIter + 1)

  Y <- Y_raw[InitMarker, ]
  outY = deconfounding(Y,K)
  Prop0 = t(outY$C$Matrix)
  allProp[[1]] <- Prop0

  out_all <- csSAM::csfit(Prop0, t(Y_raw))
  prof <- t(out_all$ghat)
  tmpmat <- prof %*% t(Prop0)
  allRMSE[1] <- sqrt(mean((t(Y_raw) - t(tmpmat)) ^ 2))

  message("+========================================+")
  message("+======= Total iterations = ",
          TotalIter, " ==========+")

  for (i in seq_len(TotalIter)) {
    message("Current iter = ", i)

    updatedInx <- DEVarSelect(Y_raw, Prop0, nMarker, bound_negative)
    Y <- Y_raw[updatedInx, ]
    Prop0 <- t(deconfounding(Y, K)$C$Matrix)

    ## avoid error by using csfit
    idx = apply(Prop0,2,function(x) sum(x) == 0)
    Prop0[,idx] = matrix(runif(ncol(Y_raw)*sum(idx),0.0001,0.0002),ncol(Y_raw),sum(idx))

    allProp[[i + 1]] <- Prop0
    out_all <- csSAM::csfit(Prop0, t(Y_raw))
    prof <- t(out_all$ghat)
    tmpmat <- prof %*% t(Prop0)
    allRMSE[i + 1] <- sqrt(mean((t(Y_raw) - t(tmpmat)) ^ 2))
  }

  min_idx <- which.min(allRMSE)
  Prop0 <- allProp[[min_idx]]

  return(list(allRMSE = allRMSE,
              allProp = allProp,
              estProp = Prop0,
              updatedInx = updatedInx
  ))
}
