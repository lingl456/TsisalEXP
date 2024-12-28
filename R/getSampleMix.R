#' Title getSampleMix
#'
#'#' Generate Mixed Samples Based on Pure Tissue References
#'
#' This function generates mixed samples by combining different tissue references with specified proportions.
#' The resulting mixed samples are simulated by first generating random proportions for each sample using
#' the Dirichlet distribution, then combining pure tissue reference panels with noise. This function can be
#' used to simulate gene expression or methylation data where the true proportions of tissue types are
#' mixed with some noise.
#'
#' @param N_sample An integer specifying the number of mixed samples to generate.
#' @param logpure_base A numeric matrix of size `p x K` representing the base values for each pure tissue.
#' The rows correspond to features (genes), and the columns correspond to tissues.
#' @param logpure_sd A numeric matrix of size `p x K` representing the standard deviation for each pure tissue.
#' @param noise_sd A numeric value specifying the standard deviation of noise to be added to the generated samples.
#' The default value is 0.1.
#'
#' @return A list containing:
#' \item{obs.Y}{A numeric matrix of size `p x N_sample` representing the generated mixed samples, where each row is a feature and each column is a sample.}
#' \item{trueProp}{A numeric matrix of size `N_sample x K` representing the true proportions of each tissue in each sample.}
#' \item{tmp}{A matrix containing the individual tissue contributions for each sample before noise is added.}
#'
#' @export
#'
#' @examples # Example of generating mixed samples
#' N_sample = 10
#' logpure_base = matrix(runif(100), nrow = 10, ncol = 10)  # Random pure tissue base
#' logpure_sd = matrix(runif(100) * 0.1, nrow = 10, ncol = 10)  # Random standard deviation for pure tissue
#' result = getSampleMix(N_sample, logpure_base, logpure_sd)
#' head(result$obs.Y)  # View the generated mixed samples
getSampleMix <- function(N_sample, logpure_base, logpure_sd, noise_sd = 0.1){
  K = ncol(logpure_base)
  p = nrow(logpure_base)

  ## get proportions:
  trueProp = getProportion(N_sample, L = K)

  ## get mix:
  obs.Y = matrix(0, p, N_sample)
  alltmp = matrix(0, p, N_sample*K)

  for(n in 1:N_sample){
    tmp = getOnePureRefPanel(logpure_base, logpure_sd)
    obs.Y[,n] = tmp %*% trueProp[n,] + rnorm(p, 0, noise_sd)
    alltmp[,seq(n, n+N_sample*(K-1), by=N_sample)] = tmp
  }

  obs.Y[obs.Y<0] <- 0.01
  rownames(obs.Y) = rownames(logpure_base)
  return(list(obs.Y = obs.Y, trueProp = trueProp,tmp = alltmp))
}
