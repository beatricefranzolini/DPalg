################################################################################
## this code runs the numerical experiment with perturbed zipf
################################################################################
rm(list = ls())

library(rstudioapi) # version 0.17.1
#set working directory to Source file directory:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("DPalg.R")

sample_sizes = c(150, 300, 600, 1500, 3000)
num_simul    = length(sample_sizes)

################################################################################
## Sample from perturbed zipf and run chains for all input sizes
################################################################################

data_zipf    = vector(mode = "list", length = num_simul)
c_true_zipf  = vector(mode = "list", length = num_simul)

SS_zipf           = vector(mode = "list", length = num_simul)
BGS_10_zipf       = vector(mode = "list", length = num_simul) 

set.seed(0)
i_count = 0
Tot     = 10000

for(n in sample_sizes){
  i_count = i_count + 1
  
  c_true_zipf[[i_count]] = sample(500, size = n, replace = TRUE, 
                                  prob = (1:500)^(-2))
  
  Y = rnorm(n, mean = c(c_true_zipf[[i_count]]*3, sd = 1))
  Y = Y - mean(Y)
  
  data_zipf[[i_count]] = Y 
  
  c_init = kmeans(Y, centers = 5)$cluster
  
  SS_zipf[[i_count]] = 
    dp_slice_mixmodel_normal_normal(Y, seed = 0, c_init = c_init, Tot = Tot)
  message(paste("ss for n =", n, "completed in", SS_zipf[[i_count]]$time))
  
  BGS_10_zipf[[i_count]] = 
    dp_BGS_mixmodel_normal_normal(Y, seed = 0, c_init = c_init, Tot = Tot)
  message(paste("BGS_10 for n =", n, "completed in", BGS_10_zipf[[i_count]]$time))
  
}

################################################################################
## Evaluate Log likelihoods 
################################################################################

i_count = 0

for(n in sample_sizes){
  
  i_count = i_count + 1
  
  SS_zipf[[i_count]]$loglik_trace = 
    dp_loglik_trace_mixmodel_normal_normal(data_zipf[[i_count]], SS_zipf[[i_count]])
  
  BGS_10_zipf[[i_count]]$loglik_trace = 
    dp_loglik_trace_mixmodel_normal_normal(data_zipf[[i_count]], BGS_10_zipf[[i_count]])
  
  message(paste("Loglik evalutations for", n, "sized dataset type 2 completed"))
} 

################################################################################
## Compute effective Sample size per second, time per 1000 iter in seconds and
#  rand indexes of the posterior point estimates
################################################################################

set.seed(0)
burn_inplus1 = floor(Tot / 2) + 1

ESS_per_second       = matrix(NA, nrow = 2, ncol = num_simul)
ESS_per_second_H     = matrix(NA, nrow = 2, ncol = num_simul)
Time_per_1000_in_sec = matrix(NA, nrow = 2, ncol = num_simul)
RI                   = matrix(NA, nrow = 2, ncol = num_simul)

for(i_count in 1:num_simul){
  Time_per_1000_in_sec[1,i_count] = 
    as.numeric(SS_zipf[[i_count]]$time, units = "secs") / (Tot / 1000)
  
  ESS_per_second[1,i_count] = ESS(SS_zipf[[i_count]]$loglik_trace[burn_inplus1:Tot]) / 
    Time_per_1000_in_sec[1,i_count] / ( (Tot - floor(Tot / 2)) / 1000)
  
  ESS_per_second_H[1,i_count] = ESS(SS_zipf[[i_count]]$H[burn_inplus1:Tot]) / 
    Time_per_1000_in_sec[1,i_count] / ( (Tot - floor(Tot / 2)) / 1000)
  
  temp = salso(do.call(rbind,SS_zipf[[i_count]]$c_samples[seq((burn_inplus1),Tot, by = 10)]))
  RI[1,i_count] = rand.index(temp, c_true_zipf[[i_count]])
}


for(i_count in 1:num_simul){
  Time_per_1000_in_sec[2,i_count] = 
    as.numeric(BGS_10_zipf[[i_count]]$time, units = "secs") / (Tot / 1000)
  
  ESS_per_second[2,i_count] = ESS(BGS_10_zipf[[i_count]]$loglik_trace[burn_inplus1:Tot]) / 
    Time_per_1000_in_sec[2,i_count] / ( (Tot - floor(Tot / 2)) / 1000)
  
  ESS_per_second_H[2,i_count] = ESS(BGS_10_zipf[[i_count]]$H[burn_inplus1:Tot]) / 
    Time_per_1000_in_sec[2,i_count] / ( (Tot - floor(Tot / 2)) / 1000)
  
  temp = salso(do.call(rbind,BGS_10_zipf[[i_count]]$c_samples[seq((burn_inplus1),Tot, by = 10)]))
  RI[2,i_count] = rand.index(temp, c_true_zipf[[i_count]])
}


Time_per_1000_in_sec
ESS_per_second
ESS_per_second_H
RI

################################################################################
## Compute rand-index per iteration and plot box-plot
################################################################################

set.seed(0)
i_count = 0

for(n in sample_sizes){
  
  i_count = i_count + 1
  c_true = c_true_zipf[[i_count]]
  
  SS_zipf[[i_count]]$ri = 
    unlist(lapply(SS_zipf[[i_count]]$c_samples[5001:10000], 
                  function(x) rand.index(c_true,x)))
  
  BGS_10_zipf[[i_count]]$ri = 
    unlist(lapply(BGS_10_zipf[[i_count]]$c_samples[5001:10000], 
                  function(x) rand.index(c_true,x)))
  
  message(paste("Rand Index evalutations for", n, "sized dataset completed"))
} 

# prepare data to plot:
SS_zipf_RI_150  = SS_zipf[[1]]$ri
SS_zipf_RI_300  = SS_zipf[[2]]$ri 
SS_zipf_RI_600  = SS_zipf[[3]]$ri 
SS_zipf_RI_1500 = SS_zipf[[4]]$ri 
SS_zipf_RI_3000 = SS_zipf[[5]]$ri 

BGS_10_zipf_RI_150  = BGS_10_zipf[[1]]$ri
BGS_10_zipf_RI_300  = BGS_10_zipf[[2]]$ri 
BGS_10_zipf_RI_600  = BGS_10_zipf[[3]]$ri 
BGS_10_zipf_RI_1500 = BGS_10_zipf[[4]]$ri 
BGS_10_zipf_RI_3000 = BGS_10_zipf[[5]]$ri 

prefixes = c("SS", "BGS_10")
labels   = c("Slice", "BGS") # labels for the legend

# build a data frame by looping over sizes and prefixes
df_list = list()
for(j in seq_along(sample_sizes)){
  N = sample_sizes[j]
  df_list[[j]] = do.call(rbind, lapply(seq_along(prefixes), function(i) {
    pref <- prefixes[i]
    varname <- paste0(pref, "_zipf_RI_", N)
    if (!exists(varname)) {
      stop("Variable not found: ", varname)
    }
    data.frame(
      sample_size = N,
      algorithm   = labels[i],
      value       = get(varname)
    )
  }))
}
# combine into one data.frame
df <- do.call(rbind, df_list)

# ensure factor ordering
df$sample_size = factor(df$sample_size, levels = sample_sizes)
df$algorithm   = factor(df$algorithm,   levels = labels)

# plot
ggplot(df, aes(x = sample_size, y = value, fill = algorithm)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    linewidth = 0.4,
    colour = "black",
    alpha = 0.85,
    outlier.size = 0.5
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    position = position_dodge(width = 0.8),
    size = 2,
    shape = 21,
    colour = "black"
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x    = "Sample Size (n)",
    y    = "Rand Index",
    fill = "Algorithm"
  ) + coord_cartesian(ylim = c(0.5,1)) +
  theme_bw(base_size = 14)

################################################################################
## plot posterior co-clustering matrices
################################################################################

i_count = 1
temp = psm(do.call(rbind,SS_zipf[[i_count]]$c_samples[5001:10000]))
image(temp, axes=F)
temp = psm(do.call(rbind,BGS_10_zipf[[i_count]]$c_samples[5001:10000]))
image(temp, axes=F)
temp = psm(c_true_zipf[[i_count]])
image(temp, axes=F)
i_count = 2
temp = psm(do.call(rbind,SS_zipf[[i_count]]$c_samples[5001:10000]))
image(temp, axes=F)
temp = psm(do.call(rbind,BGS_10_zipf[[i_count]]$c_samples[5001:10000]))
image(temp, axes=F)
temp = psm(c_true_zipf[[i_count]])
image(temp, axes=F)
i_count = 3
temp = psm(do.call(rbind,SS_zipf[[i_count]]$c_samples[5001:10000]))
image(temp, axes=F)
temp = psm(do.call(rbind,BGS_10_zipf[[i_count]]$c_samples[5001:10000]))
image(temp, axes=F)
temp = psm(c_true_zipf[[i_count]])
image(temp, axes=F)

