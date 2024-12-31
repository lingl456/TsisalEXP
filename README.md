# Analyses of gene expression data from heterogeneous samples with TsisalEXP

TsisalEXP is an R package designed for the analyses of gene expression data from complex, heterogeneous tissues.  It is designed for the analyses of gene expression data from heterogeneous tissues, which is a mixture of different cell types.

TsisalEXP offers functions for detecting cell-type specific differential expression or Estimating the number of cell types, and improving reference-free deconvolution based on cross-cell type differential analysis.   TsisalEXP implements to reference-free deconvolution with prior knowledge.

In this readme file, we briefly present how to install TsisalEXP package through GitHub. For detailed usage of TsisalEXP, please refer to the vignette file.


# Installation and quick start

## Install TsisalEXP

To install this package, start R (version "3.6") and enter:

```
if (!require("devtools", quietly = TRUE))
      install.packages("devtools")
    
   install.packages("TsisalEXP")

   library(TsisalEXP)
```

## How to get help for TsisalEXP

Any TsisalEXP questions should be posted to the GitHub Issue section of TsisalEXP homepage at [Link](https://github.com/lingl456/TsisalEXP/issues).

## quick start

In the TsisalEXP package, we provide a toy simulated dataset that you can directly generate using the following code:

```
sim_data <- getSampleMix(N_sample = Nsample, logpure, logpure.sd, noise_sd = 0.1)
  Y.raw <- sim_data$obs.Y
  trueProp = sim_data$trueProp
  trueExp = sim_data$tmp
  K = ncol(trueProp)
```
performs a deconvolution process to estimate the proportions of components in a mixture model and evaluates the quality of these estimations

```
out.TsisalExp <- csDeconv(Y.raw, K = 4, TotalIter = 5, bound_negative = TRUE)
  res.TsisalExp = mysisal(Y.raw[out.TsisalExp$updatedInx,], K = 4,topN = 50)
  estProp_TsisalExp = t(apply(res.TsisalExp$estProp,1,function(x) x/sum(x)))
  idx = apply(estProp_TsisalExp,2,function(x) sum(x) == 0)
  estProp_TsisalExp[,idx] = runif(Nsample,0.0001,0.0002)
  estProp_TsisalExp <- GetPropAligned(input = estProp_TsisalExp, reference = trueProp)
  MC_TsisalExp = mean(diag(cor(estProp_TsisalExp,trueProp)))
  MAE_TsisalExp = mean(colSums(abs(estProp_TsisalExp - trueProp)) / nrow(trueProp))
```

Some explanations about the parameters:

N_sample: Number of samples to generate.

logpure: Log-transformed pure components or distributions.

logpure.sd: Standard deviation for these pure components.

noise_sd: Standard deviation of the noise added to the observed data.

obs.Y: The observed data after adding noise.

trueProp: True proportions of each pure component in the mixture.

tmp: Additional related data, possibly the noiseless mixed data.

Y.raw: Observed noisy data.

K = 4: Number of components to estimate in the mixture.

TotalIter = 5: Number of iterations for the algorithm.

bound_negative = TRUE: A constraint to prevent negative values during the deconvolution.
