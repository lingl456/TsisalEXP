#' Title sisal
#'
#'#' SISAL: Sparse Inverse Simultaneous Component Analysis
#'
#' This function performs Sparse Inverse Simultaneous Component Analysis (SISAL) for tissue deconvolution.
#' It estimates tissue-specific gene expression profiles (endpoints) and the corresponding tissue proportions.
#' @param Y A numeric matrix of size L x N representing the gene expression data, where L is the number of genes and N is the number of samples.
#' @param p The number of tissue types/components to estimate. This should be a positive integer less than or equal to L (the number of genes).
#' @param iters An integer specifying the number of iterations for the optimization process (default is 80).
#' @param tau A regularization parameter controlling the sparsity of the solution (default is 1).
#' @param mu A regularization parameter influencing the stability of the solution (default is p * 1000 / ncol(Y)).
#' @param spherize A logical value indicating whether to apply spherization to the data (default is FALSE).
#' @param tol A numeric value specifying the convergence tolerance (default is 1e-2).
#' @param m0 An optional initial matrix for tissue proportions. If `NULL`, the matrix is computed using the VCA method (default is NULL).
#' @param verbose A logical value indicating whether to print progress messages during the iterations (default is FALSE).
#' @param returnPlot A logical value indicating whether to return a plot of the iteration steps (default is TRUE).
#' @param nonNeg A logical value indicating whether to enforce non-negativity constraints on the solution (default is FALSE).
#'
#' @return A list containing:
#' \describe{
#'   \item{endpoints}{A matrix of estimated tissue-specific components (endpoints), with columns representing the tissue types and rows representing the genes.}
#'   \item{endpointsProjection}{A matrix of projections of the data onto the tissue-specific components.}
#'   \item{projection}{A matrix of the estimated tissue proportions (projection matrix).}
#'   \item{shift}{A vector representing the mean-centered data.}
#'   \item{singValues}{A vector of singular values from the singular value decomposition (SVD).}
#'   \item{distances}{A vector of distances from each sample to the estimated tissue components (endpoints).}
#'   \item{reduced}{The reduced data matrix after dimensionality reduction.}
#'   \item{q}{The final transformation matrix used in the optimization process.}
#'   \item{plotObj}{Optional: A data frame containing the plot of iteration steps, if `returnPlot` is TRUE.}
#' }
#' @export
#'
#' @examples # Example usage of the sisal function
#' data <- matrix(rnorm(1000), nrow = 100, ncol = 10) # Random example data (100 genes x 10 samples)
#' result <- sisal(Y = data, p = 3, iters = 50, tau = 1, mu = 300)
sisal <- function(Y, p, iters = 80, tau = 1,mu = p * 1000 / ncol(Y),spherize = F,
                  tol = 1e-2, m0 = NULL, verbose=F,returnPlot = T, nonNeg = F){
  rnames <- rownames(Y)
  L <- nrow(Y)
  N <- ncol(Y)

  if (L < p) stop("Insufficient number of columns in y")

  # local stuff

  slack <- 1e-3
  energyDecreasing <- 0
  fValBack <- Inf
  lamSphe <- 1e-8
  lamQuad <- 1e-6
  ALiters <- 4
  flaged <- 0

  hinge <- function(Y) {
    return(pmax(-Y, 0))
  }

  softNeg <- function(Y, tau) {
    z <- pmax(abs(Y + tau / 2) - tau / 2, 0)
    z <- z / (z + tau / 2) * (Y + tau / 2)
    return(z)
  }

  ## At first we are getting the affine set
  Ymean <- apply(Y, 1, mean)
  ym <- matrix(Ymean, ncol=1)
  Y <- Y - Ymean
  svdObj <- fast.svd(Y)
  Up <- svdObj$u[, 1:(p-1)]
  proj <- Up %*% t(Up)
  D <- svdObj$d[1:(p-1)]

  Y <- proj %*% Y

  Y <- Y + Ymean
  YmeanOrtho <- ym - proj %*% ym
  Up <- cbind(Up, YmeanOrtho / (sqrt(sum(YmeanOrtho ^ 2))))
  singValues <- D

  Y <- t(Up) %*% Y

  ## spherizing
  if (spherize) {
    Y <- Up %*% Y
    Y <- Y - Ymean
    C <- diag(1 / sqrt(D + lamSphe))
    IC <- solve(C)
    Y <- C %*% t(Up[, 1:(p-1)]) %*% Y
    Y <- rbind(Y, 1)
    Y <- Y / sqrt(p)
  }

  ## Init

  if (is.null(m0)) {
    Mvca <- vca(Y, p, verbose=verbose)
    M <- Mvca
    Ym <- apply(M, 1, mean)
    dQ <- M - Ym
    M <- M + p * dQ
  } else {
    M <- m0
    M <- M - Ymean
    M <- Up[, 1:(p-1)] %*% t(Up[, 1:(p - 1)]) %*% M
    M <- M + Ymean
    M <- t(Up) %*% M
    if (spherize) {
      M <- Up %*% M - Ymean
      M <- C %*% t(Up[, 1:(p - 1)]) %*% M
      M[p, ] <- 1
      M <- M / sqrt(p)
    }
  }

  if (returnPlot) {
    toPlot <- data.frame(
      x = Y[1, ],
      y = Y[2, ],
      type = "data point",
      iter = NA,
      tau = NA
    )
    starting <- data.frame(
      x = M[1, ],
      y = M[2, ],
      type = "sisal",
      iter = 0,
      tau = tau
    )
    toPlot <- rbind(toPlot, starting)
  }


  Q0 <- solve(M)
  Q <- Q0

  AAT <- kronecker(Y %*% t(Y), diag(nrow=p))
  B <- kronecker(diag(nrow=p), matrix(1, nrow=1, ncol=p))
  qm <- rowSums(solve(Y %*% t(Y)) %*% Y)
  qm <- matrix(qm, ncol=1)

  H <- lamQuad * diag(nrow=p^2)
  FF <- H + mu * AAT
  IFF <- solve(FF)


  G <- IFF %*% t(B) %*% solve(B %*% IFF %*% t(B))
  qmAux <- G %*% qm
  G <- IFF - G %*% B %*% IFF

  Z <- Q %*% Y
  Bk <- 0 * Z

  fmin = Inf
  Qmin = NULL

  for (k in 1:iters) {
    IQ <- solve(Q)
    g <- -t(IQ)
    dim(g) <- c(nrow(g) * ncol(g), 1)

    q0 <- Q
    dim(q0) <- c(nrow(Q) * ncol(Q), 1)
    Q0 <- Q

    baux <- H %*% q0 - g

    if (verbose) {
      if (spherize) {
        M <- IQ * sqrt(p)
        M <- M[1:(p-1), ]
        M <- Up[, 1:(p-1)] %*% IC %*% M
        M <- M + Ymean
        M <- t(Up) %*% M
      } else {
        M <- IQ
      }
      message(sprintf("Iteration %d, simplex volume = %4f", k, abs(det(M)) / factorial(nrow(M))))
    }

    if (k == iters) {
      ALiters <- 100
    }

    while (T) {
      q <- Q
      dim(q) <- c(nrow(Q) * ncol(Q), 1)

      f0val <- -log(abs(det(Q))) + tau * sum(hinge(Q %*% Y))
      f0quad <- t(q - q0) %*% g + 0.5 * t(q - q0) %*% H %*% (q - q0) + tau * sum(hinge(Q %*% Y))

      for (i in 2:ALiters) {
        dqAux <- Z + Bk
        dtzB <- dqAux %*% t(Y)
        dim(dtzB) <- c(nrow(dtzB) * ncol(dtzB), 1)
        b <- baux + mu * dtzB
        q <- G %*% b + qmAux
        Q <- matrix(q, nrow=p, ncol=p)

        Z <- softNeg(Q %*% Y - Bk, tau / mu)

        Bk <- Bk - (Q %*% Y - Z)

      }

      fquad_tmp <- t(q - q0) %*% g + 0.5 * t(q - q0) %*% H %*% (q - q0) + tau * sum(hinge(Q %*% Y))
      fval_tmp <- -log(abs(det(Q))) + tau * sum(hinge(Q %*% Y))

      fquad <- t(q - q0) %*% g + 0.5 * t(q - q0) %*% H %*% (q - q0) + tau * sum(hinge(Q %*% Y))
      fval <- -log(abs(det(Q))) + tau * sum(hinge(Q %*% Y))

      if (f0quad >= fquad) {
        while (f0val - fval < 0) {
          Q <- (Q + Q0) / 2
          fval <- -log(abs(det(Q))) + tau * sum(hinge(Q %*% Y))
        }

        break
      }

    }

    if (returnPlot) {
      M <- solve(Q)
      toAdd <- data.frame(
        x = M[1, ],
        y = M[2, ],
        type = "sisal",
        iter = k,
        tau = tau
      )
      toPlot <- rbind(toPlot, toAdd)
      toPlot <- tibble::as_tibble(toPlot)
    }

  }

  if (nonNeg) {
    qorig <- Up %*% solve(Q)

    if (any(qorig < 0)) {
      diff <- Ymean - qorig

      shift <- sapply(1:p, function(i) {
        tmp <- qorig[, i] / diff[, i]
        ifelse(any(qorig[, i] < 0), min(tmp[qorig[, i] < 0]), 0)
      })

      shift <- -shift + (0.00001)
      qorig <- qorig + diff %*% diag(shift)
      Q <- solve(t(Up) %*% qorig)
      q <- Q
      dim(q) <- c(nrow(Q) * ncol(Q), 1)
    }
  }

  distanceToEndpoints <- apply(M, 2, function(x) {
    shifted <- Y - x
    dds <- sqrt(colSums(shifted^2))
    return(dds)
  })

  endpointsProjection <- M

  if (spherize) {
    M <- solve(Q)
    M <- M * sqrt(p)
    M <- M[1:(p - 1), ]
    M <- Up[, 1:(p - 1)] %*% IC %*% M
    M <- M + Ymean
  } else {
    M <- Up %*% solve(Q)
  }

  colnames(M) <- paste0("Pure gene ", 1:p)
  rownames(M) <- rnames

  retVal <- list(
    endpoints = M,
    endpointsProjection = endpointsProjection,
    projection = Up,
    shift = Ymean,
    singValues = singValues,
    distances = distanceToEndpoints,
    reduced = Y,
    q = Q
  )
  if (returnPlot) retVal$plotObj <- toPlot
  return(retVal)
}
