#' Title vca
#'
#' Variance Component Analysis (VCA)
#'
#' This function performs Variance Component Analysis (VCA) to extract relevant components from the input data matrix.
#' VCA is often used for decomposing variance in gene expression data to identify principal components or sources of variation.
#'
#' @param R A numeric matrix of size L x N, where L is the number of genes (rows) and N is the number of samples (columns).
#' Each element represents the expression level of a gene in a sample.
#' @param p An integer specifying the number of components to extract.
#' @param SNR Optional. A numeric value for the Signal-to-Noise Ratio. If provided, the analysis will adapt based on the SNR threshold.
#' @param verbose A logical value indicating whether to print progress messages. Default is `FALSE`.
#'
#' @return A numeric matrix of size L x p, containing the extracted components from the input data matrix R.
#' The columns of the returned matrix represent the variance components.
#' @export
#'
#' @examples # Example usage of the vca function
#' R <- matrix(rnorm(1000), nrow = 100, ncol = 10)  # Random example data (100 genes x 10 samples)
#' p <- 3  # Number of components to extract
#' vca_results <- vca(R, p)
vca <- function(R, p, SNR=NULL, verbose=F) {
  L <- nrow(R)
  N <- ncol(R)
  if (p < 0 || p > L) {
    stop("p is out of range (negative or too big)")
  }
  SNRth <- 15 + 10 * log10(p)

  if (!is.null(SNR) && (SNR < SNRth)) {
    if (verbose) message("Select the projective projection")
    d <- p - 1

    if (exists("xp")) {
      Ud <- Ud[, 1:d]
    } else {
      rm <- apply(R, 1, mean)
      R0 <- R - rm
      svdObj <- svd(R0 %*% t(R0), nu = p, nv = p)
      Ud <- svdObj$u
      xp <- t(Ud) %*% R0
    }

    Rp <- Ud %*% xp[1:d, ] + rm
    x <- xp[1:d, ]
    c <- sqrt(max(colSums(x^2)))
    y <- rbind(x, rep(c, N))
  } else {
    if (verbose) message("Select projection to p-1")
    d <- p
    Ud <- svd(R %*% t(R) / N, nu = d, nv = d)$u

    xp <- t(Ud) %*% R
    Rp <- Ud %*% xp[1:d, ]

    x <- xp
    u <- rowMeans(x)
    y <- x / matrix(kronecker(colSums(x * u), rep(1, d)), nrow=d)
  }

  ## VCA itself

  indice <- rep(0, p)
  A <- matrix(0, nrow=p, ncol=p)
  A[p, 1] <- 1

  for (i in 1:p) {
    w <- matrix(runif(p), ncol=1)
    f <- w - A %*% pseudoinverse(A) %*% w;
    f <- f / sqrt(sum(f^2))

    v <- t(f) %*% y
    indice[i] <- which.max(abs(v))
    A[, i] <- y[, indice[i]]
  }
  Ae = Rp[, indice]
  return(Ae)
}

