rm(list = ls())

library(rstudioapi) # version 0.17.1
#set working directory to Source file directory:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("DPalg.R")

sample_sizes = c(150, 300, 600, 1500, 3000, 7500, 12000)
num_simul    = length(sample_sizes)

#sample from 3 well separated clusters 
data    = vector(mode = "list", length = num_simul)
c_true  = vector(mode = "list", length = num_simul)

SS           = vector(mode = "list", length = num_simul)
SSnoAtoms    = vector(mode = "list", length = num_simul)
BGS_10       = vector(mode = "list", length = num_simul) 
BGS_n        = vector(mode = "list", length = num_simul)
CRPwithAtoms = vector(mode = "list", length = num_simul)
CRPnoAtoms   = vector(mode = "list", length = num_simul) 

set.seed(0)
i_count = 0

for(n in sample_sizes){
  i_count = i_count + 1
  
  c_true[[i_count]] = c(rep(1, floor(n/3)), 
                        rep(2, floor(n/3)), 
                        rep(3,n - 2*floor(n/3)))
  
  Y = rnorm(n, mean = c(rep(-3, floor(n/3)), 
                         rep(0, floor(n/3)), 
                         rep(3,n - 2*floor(n/3))), sd = 1)
  
  data[[i_count]] = Y
  
  c_init = kmeans(Y, centers = 5)$cluster
  
  SS[[i_count]] = 
    dp_slice_mixmodel_normal_normal(Y, seed = 0, c_init = c_init, Tot = 5000)
  message(paste("ss for n =", n, "completed in", SS[[i_count]]$time))
  
  SSnoAtoms[[i_count]] = 
    dp_slicenoAtoms_mixmodel_normal_normal(Y, seed = 0, c_init = c_init, Tot = 5000)
  message(paste("ss no atoms for n =", n, "completed in", SSnoAtoms[[i_count]]$time))
  
  BGS_10[[i_count]] = 
    dp_BGS_mixmodel_normal_normal(Y, seed = 0, c_init = c_init, Tot = 5000)
  message(paste("BGS_10 for n =", n, "completed in", BGS_10[[i_count]]$time))
  
  BGS_n[[i_count]] = 
    dp_BGS_mixmodel_normal_normal(Y, L = n, seed = 0, c_init = c_init, Tot = 5000)
  message(paste("BGS_n for n =", n, "completed in", BGS_n[[i_count]]$time))
  
  CRPwithAtoms[[i_count]] = 
    dp_CRPwithAtoms_mixmodel_normal_normal(Y, seed = 0, c_init = c_init, Tot = 5000)
  message(paste("CRPwithAtoms for n =", n, "completed in", CRPwithAtoms[[i_count]]$time))
  
  CRPnoAtoms[[i_count]] = 
    dp_CRPnoAtoms_mixmodel_normal_normal(Y, seed = 0, c_init = c_init, Tot = 5000)
  message(paste("CRPnoAtoms for n =", n, "completed in", CRPnoAtoms[[i_count]]$time))
  
}

#sample from perturbed zipf 
data_zipf    = vector(mode = "list", length = num_simul)
c_true_zipf  = vector(mode = "list", length = num_simul)
  
SS_zipf           = vector(mode = "list", length = num_simul)
CRPnoAtoms_zipf   = vector(mode = "list", length = num_simul) 
CRPwithAtoms_zipf = vector(mode = "list", length = num_simul)
BGS_10_zipf       = vector(mode = "list", length = num_simul) 
BGS_n_zipf        = vector(mode = "list", length = num_simul)

set.seed(0)
i_count = 0

for(n in sample_sizes){
  i_count = i_count + 1
  
  c_true_zipf[[i_count]] = sample(500, size = n, replace = TRUE, 
                                  prob = (1:500)^(-2))
  
  Y = rnorm(n, mean = c(c_true_zipf[[i_count]]*3, sd = 1))
  
  data_zipf[[i_count]] = Y
  
  c_init = kmeans(Y, centers = 5)$cluster
  
  SS_zipf[[i_count]] = 
    dp_slice_mixmodel_normal_normal(Y, seed = 0, c_init = c_init, Tot = 5000)
  message(paste("ss for n =", n, "completed in", SS[[i_count]]$time))
  
  CRPnoAtoms_zipf[[i_count]] = 
    dp_CRPnoAtoms_mixmodel_normal_normal(Y, seed = 0, c_init = c_init, Tot = 5000)
  message(paste("CRPnoAtoms for n =", n, "completed in", CRPnoAtoms[[i_count]]$time))
  
  CRPwithAtoms_zipf[[i_count]] = 
    dp_CRPwithAtoms_mixmodel_normal_normal(Y, seed = 0, c_init = c_init, Tot = 5000)
  message(paste("CRPwithAtoms for n =", n, "completed in", CRPwithAtoms[[i_count]]$time))
  
  BGS_10_zipf[[i_count]] = 
    dp_BGS_mixmodel_normal_normal(Y, seed = 0, c_init = c_init, Tot = 5000)
  message(paste("BGS_10 for n =", n, "completed in", BGS_10[[i_count]]$time))
  
  BGS_n_zipf[[i_count]] = 
    dp_BGS_mixmodel_normal_normal(Y, L = n, seed = 0, c_init = c_init, Tot = 5000)
  message(paste("BGS_n for n =", n, "completed in", BGS_n[[i_count]]$time))
  
}

