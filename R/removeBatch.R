#' Title removeBatch
#'
#'#' Remove Batch Effects from Data Using ComBat
#'
#' This function removes batch effects from gene expression data by applying the ComBat algorithm.
#' It uses a reference dataset (`Jaff.blood.ref`) to remove systematic differences between
#' the two data sets (`dat.raw` and `Jaff.blood.ref`).
#'
#' @param dat.raw A numeric matrix or data frame containing the raw data. The rows represent features (e.g., genes),
#'                and the columns represent samples. The samples in `dat.raw` should match the samples in `Jaff.blood.ref`.
#' @param Jaff.blood.ref A numeric matrix or data frame containing reference data to which batch effects will be removed.
#'                       The rows should correspond to the same features as in `dat.raw`. The columns 3 to 8 in `Jaff.blood.ref`
#'                       are used for the batch effect adjustment.
#'
#' @return A list with two components:
#' \item{Y.raw}{A matrix containing the gene expression data with the batch effects removed. It is the adjusted data from `dat.raw`.}
#' \item{pure}{A matrix containing the reference data with the batch effects removed. It corresponds to the adjusted data from `Jaff.blood.ref`.}
#'
#' @export
#'
#' @examples
removeBatch <- function(dat.raw, Jaff.blood.ref){

  print("Removing batch effects!")

  ol.names = intersect(rownames(Jaff.blood.ref),rownames(dat.raw))

  Ytmp = dat.raw[ol.names,]
  Reftmp = Jaff.blood.ref[ol.names,3:8]
  # NewReftmp = cbind(apply(Reftmp[,1:3],1,mean),Reftmp[,4:6])
  # Reftmp  = NewReftmp

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
