# DPalg
R codes implementing MCMC algoritms for for DP mixture of univariate Normals

Kernel:   y | mu ~ N(mu, sigma2)

Base:     mu ~ N(mu0, tau20)

DP:       G ~ DP(alpha, N(mu0, tau20)) with alpha random

**Authors**: [blinded](blinded) and [blinded](blinded)

## Citation 
Please cite the following publication if you use this repository in your research: 

[blinded](blinded) and [blinded](blinded), Complexity bounds for Dirichlet process slice samplers.

## Contents 

### Main scripts 

`_main_to_run1.R` -> runs the numerical experiment with 3 equally-sized clusters

`_main_to_run2.R` -> runs the numerical experiment with perturbed zipf

`BGS.R`           -> MCMC function to run the block Gibbs sampler

`CRPnoAtoms.R`    -> MCMC function to run the CRP with atoms

`CRPwithAtoms.R`  -> MCMC function to run the CRP with atoms

`DPalg.R`         -> this code sources the required functions and packages it is called by _main_to_run1.R and _main_to_run2.R

`Slice.R`         -> MCMC function to run the slice sampler

`SlicenoAtoms.R`  -> MCMC function to run the slice sampler with no atoms

`utils`           -> contains helper functions for the MCMCs

### R libraries needed to run the scripts
library(LaplacesDemon)  # version 16.1.6      to sample from finite Dirichlet  

library(progress)       # version 1.2.3       to draw the progress bar  

library(rstudioapi)     # version 0.17.1      to set working directory  

library(ggplot2)        # version 4.0.1       to plot  

library(patchwork)      # version 1.3.2       to plot  

library(fossil)         # version 0.4.0       to compute rand index  

library(salso)          # version 0.3.53      to compute psm and point estimate

### Questions or bugs
For bug reporting purposes, e-mail [blinded] (blinded) (author@author.me).