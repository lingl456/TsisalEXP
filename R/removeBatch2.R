#' Title removeBatch2
#'
#'#' Remove Batch Effects from TCGA Data
#'
#' This function removes batch effects from TCGA data by applying the ComBat method from the `sva` package.
#' It assumes that the raw data and reference data are matched by gene names.
#'
#' @param dat.raw A numeric matrix of raw TCGA data, with rows representing genes and columns representing samples.
#' @param pure.mean A numeric matrix of reference (pure) data, with rows representing genes and columns representing samples.
#' The reference data should have the same gene names as the raw data.
#'
#'
#' @return A list containing:
#' \describe{
#'   \item{Y.raw}{A numeric matrix of the batch-corrected TCGA data.}
#'   \item{pure}{A numeric matrix of the batch-corrected reference data.}
#' }
#'
#' @export
#'
#' @examples # Example usage:
#' result <- removeBatch2(dat.raw = tcga_data, pure.mean = reference_data)
#' batch_corrected_data <- result$Y.raw
#' corrected_reference_data <- result$pure
removeBatch2 <- function(dat.raw,pure.mean){

  print("Removing batch effects!")

  ol.names = intersect(rownames(pure.mean),rownames(dat.raw))

  Ytmp = dat.raw[ol.names,]
  Reftmp = pure.mean[ol.names,]

  ## remove batch effect
  library(sva)
  batch = c(rep(1,ncol(Ytmp)),rep(2,ncol(Reftmp)))
  batch = as.numeric(batch)
  Ytotal = as.matrix(cbind(Ytmp,Reftmp))
  names(batch) = colnames(Ytotal)
  combat_edata = ComBat(dat=Ytotal, batch=batch)

  Y.raw = combat_edata[,1:ncol(Ytmp)]
  pure = combat_edata[,ncol(Ytmp)+1:ncol(Reftmp)]

  print("Finish - Removing batch effects!")

  return(list(Y.raw = Y.raw,
              pure = pure))
}

