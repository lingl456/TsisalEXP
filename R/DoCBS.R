#' Title
#'
#'#' Perform CBS (Constrained Basis Selection) Using Support Vector Machines
#'
#' This function performs CBS (Constrained Basis Selection) using Support Vector Machines (SVM) to select the optimal basis
#' for reconstructing data. The function estimates the coefficients for each sample using a linear SVM with a specified
#' \code{nu} parameter, and then selects the best \code{nu} based on the root mean squared error (RMSE) between the observed
#' and reconstructed values.
#' @param Y A numeric matrix of size N x P, where N is the number of samples and P is the number of target variables (e.g., genes).
#' @param W A numeric matrix of size N x Q, where N is the number of samples and Q is the number of basis (e.g., features).
#' @param nu.v A numeric vector of values for the \code{nu} parameter in the SVM model, default is \code{c(0.25, 0.5, 0.75)}.
#' The \code{nu} parameter controls the upper bound on the fraction of margin errors and the lower bound on the fraction of support vectors.
#'
#' @return A numeric matrix of size P x Q, where P is the number of target variables and Q is the number of basis, containing
#' the estimated coefficients for each target variable and basis after selecting the optimal \code{nu}.
#' @export
#'
#' @examples # Example usage of DoCBS
#' Y <- matrix(rnorm(100), nrow = 10, ncol = 10)  # Example target matrix (10 samples x 10 variables)
#' W <- matrix(rnorm(100), nrow = 10, ncol = 5)   # Example basis matrix (10 samples x 5 features)
#' H <- DoCBS(Y, W, nu.v = c(0.25, 0.5, 0.75))
DoCBS <- function(Y, W, nu.v = c(0.25, 0.5, 0.75)) {
  require(e1071)

  if (is.vector(W)) {
    W = matrix(W)
  }

  est.lm <- list()

  nui <- 1

  for (nu in nu.v) {
    est.m <- matrix(nrow = ncol(Y), ncol = ncol(W))

    for (s in 1:ncol(Y)) {
      svm.o <-
        svm(
          x = W,
          y = Y[, s],
          scale = TRUE,
          type = "nu-regression",
          kernel = "linear",
          nu = nu
        )

      coef.v <- t(svm.o$coefs) %*% svm.o$SV

      coef.v[which(coef.v < 0)] <- 0

      total <- sum(coef.v)

      coef.v <- coef.v / total

      est.m[s, ] <- coef.v

    }
    est.lm[[nui]] <- est.m

    nui <- nui + 1

  }

  H = matrix(nrow = ncol(Y), ncol = ncol(W))
  #### select best nu
  rmse.m <- matrix(NA, nrow = ncol(Y), ncol = length(nu.v))

  for (nui in 1:length(nu.v)) {
    reconst.m <- W %*% t(est.lm[[nui]])

    for (s in 1:ncol(Y)) {
      rmse.m[s, nui] <- sqrt(mean((Y[, s] - reconst.m[, s]) ^ 2))

    }
  }
  nu.idx <- apply(rmse.m, 1, which.min)

  H <- est.m

  for (s in 1:ncol(Y)) {
    H[s, ] <- est.lm[[nu.idx[s]]][s, ]

  }
  return(H)

}
