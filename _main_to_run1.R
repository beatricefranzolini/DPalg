################################################################################
## this code runs the numerical experiment with 3 equally-sized clusters
################################################################################
rm(list = ls())

library(rstudioapi) # version 0.17.1
#set working directory to Source file directory:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("DPalg.R")

sample_sizes = c(150, 300, 600, 1500, 3000, 7500, 12000)
num_simul    = length(sample_sizes)

################################################################################
## Sample from 3 equally-sized clusters and run chains for all input sizes
################################################################################

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
Tot     = 10000
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
    dp_slice_mixmodel_normal_normal(Y, seed = 0, c_init = c_init, Tot = Tot)
  message(paste("ss for n =", n, "completed in", SS[[i_count]]$time))
  
  SSnoAtoms[[i_count]] = 
    dp_slicenoAtoms_mixmodel_normal_normal(Y, seed = 0, c_init = c_init, Tot = Tot)
  message(paste("ss no atoms for n =", n, "completed in", SSnoAtoms[[i_count]]$time))
  
  BGS_10[[i_count]] = 
    dp_BGS_mixmodel_normal_normal(Y, seed = 0, c_init = c_init, Tot = Tot)
  message(paste("BGS_10 for n =", n, "completed in", BGS_10[[i_count]]$time))
  
  BGS_n[[i_count]] = 
    dp_BGS_mixmodel_normal_normal(Y, L = n, seed = 0, c_init = c_init, Tot = Tot)
  message(paste("BGS_n for n =", n, "completed in", BGS_n[[i_count]]$time))
  
  CRPwithAtoms[[i_count]] = 
    dp_CRPwithAtoms_mixmodel_normal_normal(Y, seed = 0, c_init = c_init, Tot = Tot)
  message(paste("CRPwithAtoms for n =", n, "completed in", CRPwithAtoms[[i_count]]$time))
  
  CRPnoAtoms[[i_count]] = 
    dp_CRPnoAtoms_mixmodel_normal_normal(Y, seed = 0, c_init = c_init, Tot = Tot)
  message(paste("CRPnoAtoms for n =", n, "completed in", CRPnoAtoms[[i_count]]$time))
  
}

################################################################################
## Evaluate Log likelihoods 
################################################################################

i_count = 0

for(n in sample_sizes){
  
  i_count = i_count + 1
  
  SS[[i_count]]$loglik_trace = 
    dp_loglik_trace_mixmodel_normal_normal(data[[i_count]], SS[[i_count]])
  
  SSnoAtoms[[i_count]]$loglik_trace = 
    dp_loglik_trace_mixmodel_normal_normal(data[[i_count]], SSnoAtoms[[i_count]])
  
  BGS_10[[i_count]]$loglik_trace = 
    dp_loglik_trace_mixmodel_normal_normal(data[[i_count]], BGS_10[[i_count]])
  
  BGS_n[[i_count]]$loglik_trace = 
    dp_loglik_trace_mixmodel_normal_normal(data[[i_count]], BGS_n[[i_count]])
  
  CRPwithAtoms[[i_count]]$loglik_trace = 
    dp_loglik_trace_mixmodel_normal_normal(data[[i_count]], CRPwithAtoms[[i_count]])
  
  CRPnoAtoms[[i_count]]$loglik_trace = 
    dp_loglik_trace_mixmodel_normal_normal(data[[i_count]], CRPnoAtoms[[i_count]])
  
  message(paste("Loglik evalutations for", n, "sized dataset type 1 completed"))
} 

################################################################################
## Compute effective Sample size per second, time per 1000 iter in seconds and
#  rand indexes of the posterior point estimates
################################################################################

set.seed(0)

ESS_per_second       = matrix(NA, nrow = 6, ncol = num_simul)
Time_per_1000_in_sec = matrix(NA, nrow = 6, ncol = num_simul)
RI                   = matrix(NA, nrow = 6, ncol = num_simul)

burn_inplus1 = floor(Tot / 2) + 1 #burnin

for(i_count in 1:num_simul){
  Time_per_1000_in_sec[1,i_count] = 
    as.numeric(SS[[i_count]]$time, units = "secs") / (Tot / 1000)
  
  ESS_per_second[1,i_count] = ESS(SS[[i_count]]$loglik_trace[burn_inplus1:Tot]) / 
    Time_per_1000_in_sec[1,i_count] / ( (Tot - floor(Tot / 2)) / 1000)
  
  temp = salso(do.call(rbind,SS[[i_count]]$c_samples[seq((burn_inplus1),Tot, by = 10)]))
  RI[1,i_count] = rand.index(temp, c_true[[i_count]])
}

for(i_count in 1:num_simul){
  if(!is.null(SSnoAtoms[[i_count]]$loglik_trace[Tot])){
    Time_per_1000_in_sec[2,i_count] = 
      as.numeric(SSnoAtoms[[i_count]]$time, units = "secs") / (Tot / 1000)
    
    ESS_per_second[2,i_count] = ESS(SSnoAtoms[[i_count]]$loglik_trace[burn_inplus1:Tot]) / 
      Time_per_1000_in_sec[2,i_count] / ( (Tot - floor(Tot / 2)) / 1000)
    
    temp = salso(do.call(rbind,SSnoAtoms[[i_count]]$c_samples[seq((burn_inplus1),Tot, by = 10)]))
    RI[2,i_count] = rand.index(temp, c_true[[i_count]])
  }
}

for(i_count in 1:num_simul){
  Time_per_1000_in_sec[3,i_count] = 
    as.numeric(BGS_10[[i_count]]$time, units = "secs") / (Tot / 1000)
  
  ESS_per_second[3,i_count] = ESS(BGS_10[[i_count]]$loglik_trace[burn_inplus1:Tot]) / 
    Time_per_1000_in_sec[3,i_count] / ( (Tot - floor(Tot / 2)) / 1000)
  
  temp = salso(do.call(rbind,BGS_10[[i_count]]$c_samples[seq((burn_inplus1),Tot, by = 10)]))
  RI[3,i_count] = rand.index(temp, c_true[[i_count]])
}

for(i_count in 1:num_simul){
  if(!is.null(BGS_n[[i_count]]$loglik_trace[Tot])){
    Time_per_1000_in_sec[4,i_count] = 
      as.numeric(BGS_n[[i_count]]$time, units = "secs") / (Tot / 1000)
    
    ESS_per_second[4,i_count] = ESS(BGS_n[[i_count]]$loglik_trace[burn_inplus1:Tot]) / 
      Time_per_1000_in_sec[4,i_count] / ( (Tot - floor(Tot / 2)) / 1000)
    
    temp = salso(do.call(rbind,BGS_n[[i_count]]$c_samples[seq((burn_inplus1),Tot, by = 10)]))
    RI[4,i_count] = rand.index(temp, c_true[[i_count]])
  }
}

for(i_count in 1:num_simul){
  if(!is.null(CRPwithAtoms[[i_count]]$loglik_trace[Tot])){
    Time_per_1000_in_sec[5,i_count] = 
      as.numeric(CRPwithAtoms[[i_count]]$time, units = "secs") / (Tot / 1000)
    
    ESS_per_second[5,i_count] = ESS(CRPwithAtoms[[i_count]]$loglik_trace[burn_inplus1:Tot]) / 
      Time_per_1000_in_sec[5,i_count] / ( (Tot - floor(Tot / 2)) / 1000)
    
    temp = salso(do.call(rbind,CRPwithAtoms[[i_count]]$c_samples[seq((burn_inplus1),Tot, by = 10)]))
    RI[5,i_count] = rand.index(temp, c_true[[i_count]])
  }
}

for(i_count in 1:num_simul){
  if(!is.null(CRPnoAtoms[[i_count]]$loglik_trace[Tot])){
    Time_per_1000_in_sec[6,i_count] = 
      as.numeric(CRPnoAtoms[[i_count]]$time, units = "secs") / (Tot / 1000)
    
    ESS_per_second[6,i_count] = ESS(CRPnoAtoms[[i_count]]$loglik_trace[burn_inplus1:Tot]) / 
      Time_per_1000_in_sec[6,i_count] / ( (Tot - floor(Tot / 2)) / 1000)
    
    temp = salso(do.call(rbind,CRPnoAtoms[[i_count]]$c_samples[seq((burn_inplus1),Tot, by = 10)]))
    RI[6,i_count] = rand.index(temp, c_true[[i_count]])
  }
}

#print tables
Time_per_1000_in_sec
ESS_per_second
RI

################################################################################
## Plot posterior similarity matrix for $n=600$
################################################################################

i_count = 3
temp = psm(do.call(rbind,SS[[i_count]]$c_samples[5001:10000]))
image(temp, axes=F)
temp = psm(do.call(rbind,SSnoAtoms[[i_count]]$c_samples[5001:10000]))
image(temp, axes=F)
temp = psm(do.call(rbind,BGS_10[[i_count]]$c_samples[5001:10000]))
image(temp, axes=F)
temp = psm(do.call(rbind,BGS_n[[i_count]]$c_samples[5001:10000]))
image(temp, axes=F)
temp = psm(do.call(rbind,CRPwithAtoms[[i_count]]$c_samples[5001:10000]))
image(temp, axes=F)
temp = psm(do.call(rbind,CRPnoAtoms[[i_count]]$c_samples[5001:10000]))
image(temp, axes=F)
temp = psm(c_true[[i_count]])
image(temp, axes=F)