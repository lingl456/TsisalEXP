#' Title getCellNumber
#'
#'#' Estimate the Optimal Number of Cells Using AIC
#'
#' This function estimates the optimal number of cell types by applying a deconvolution method to a given gene expression matrix
#' and comparing the resulting models using the Akaike Information Criterion (AIC). The function iterates over a range of possible
#' cell numbers and selects the one with the lowest AIC value.
#' @param Y.raw A numeric matrix or a SummarizedExperiment object representing gene expression data.
#' The rows represent genes (features) and columns represent samples.
#' @param possibleCellNumber A numeric vector specifying the possible number of cell types to evaluate, default is 2:15.
#'
#' @return A list with the following elements:
#' \item{bestK}{The optimal number of cell types based on the lowest AIC value.}
#' \item{allAIC}{A numeric vector of AIC values for each tested cell number.}
#' @export
#'
#' @examples # Example usage of getCellNumber
#' Y_raw <- matrix(rnorm(1000), nrow = 50, ncol = 20)  # Example gene expression matrix (50 genes x 20 samples)
#' result <- getCellNumber(Y.raw = Y_raw, possibleCellNumber = 2:10)
#' print(result$bestK)
#' print(result$allAIC)
getCellNumber <- function(Y.raw, possibleCellNumber = 2:15){

  allAIC = c()
  for(K in possibleCellNumber){
    Y_raw = Y.raw
    K = 10
    nMarker = 1000
    TotalIter = 5
    InitMarker = NULL
    bound_negative = TRUE


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
    estProp = Prop0
    estProp = mysisal(Y.raw[updatedInx,], K = K, topN = 50)$estProp
    aic = compute_aic(estProp,Y.raw[updatedInx,])
    allAIC = append(allAIC,aic)
  }
  bestK = possibleCellNumber[which.min(allAIC)]
  return(list(bestK = bestK,
              allAIC = allAIC))
}

