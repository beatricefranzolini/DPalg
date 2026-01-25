#this code source the required functions and packages 
library(LaplacesDemon)  # version 16.1.6         to sample from finite Dirichlet
library(progress)       # version 1.2.3          to draw the progress bar
#library(rstudioapi)    # version 0.17.1 (up in _main_to_run)

source("Slice.R")           #posterior Slice Sampler MCMC 
source("SlicenoAtoms.R")    #posterior Slice Sampler MCMC with marginalized atoms
source("CRPwithAtoms.R")    #posterior marginal MCMC with atom sampling 
source("CRPnoAtoms.R")      #posterior marginal MCMC with marginalized atoms 
source("BGS.R")             #Block Gibbs sampler
source("utils.R")           #Likelihoods eval & atoms sampling 
