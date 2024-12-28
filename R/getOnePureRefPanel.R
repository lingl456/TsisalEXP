#' Title getOnePureRefPanel
#'
#'#' Generate Synthetic Reference Panel for Pure Tissues
#'
#' This function generates a synthetic reference panel for gene expression or DNA methylation data.
#' It uses a normal distribution with the log-transformed means (`logpure_base`) and standard deviations
#' (`logpure_sd`) for each feature (e.g., gene or CpG site) in each tissue type to simulate the pure tissue values.
#' The generated values are then exponentiated to return values in the original scale (i.e., non-log-transformed).
#'
#' @param logpure_base A numeric matrix of log-transformed mean values for each feature (rows) across pure tissue types (columns).
#' @param logpure_sd A numeric matrix of log-transformed mean values for each feature (rows) across pure tissue types (columns).
#'
#' @return A numeric matrix of synthetic reference panel values, where each column represents one pure tissue and
#' each row corresponds to a feature (e.g., gene or CpG site). The values are in the original scale (exponentiated from log-transformed).
#'
#' @export
#'
#' @examples #' # Example input matrices for log-transformed mean and standard deviation
#' logpure_base <- matrix(runif(1000), nrow = 100, ncol = 10)  # 100 features, 10 tissues
#' logpure_sd <- matrix(runif(1000), nrow = 100, ncol = 10)  # Standard deviations
#'
#' # Generate one synthetic reference panel
#' synthetic_ref_panel <- getOnePureRefPanel(logpure_base, logpure_sd)
#'
#' # View the result
#' head(synthetic_ref_panel)
getOnePureRefPanel = function(logpure_base, logpure_sd){

  N_feature = dim(logpure_base)[1] # number of features
  L = dim(logpure_base)[2] # number of pure tissues

  tissue = matrix(0,N_feature,L)
  for(i in 1:L){
    tissue[,i] = exp(rnorm(N_feature, logpure_base[,i], abs(logpure_sd[,i])))
  }

  return(tissue)
}
