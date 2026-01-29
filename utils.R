################################################################################
## this code contain helper functions for the MCMCs
#  for DP mixture of univariate Normals
#  Kernel:   y | mu ~ N(mu, sigma2)
#  Base:     mu ~ N(mu0, tau20)
#  DP:       G ~ DP(alpha, N(mu0, tau20)) with alpha random
################################################################################

# ---- helper: compute log lik at each iteration ----
# used by _main_to_run1.R and _main_to_run2.R
dp_loglik_trace_mixmodel_normal_normal <- function(Y, fit, hyper = c(1, 0, 1)) {
  sd_y   = sqrt(hyper[1])
  
  # Make sure Y is a plain numeric vector
  if (is.list(Y)) Y = unlist(Y)
  Y = as.numeric(Y)
  
  S = length(fit$c_samples)
  if (is.null(fit$phis[[S]])){return(NULL)}
  
  ll = numeric(S)
  
  for (t in seq_len(S)) {
    c_t    = fit$c_samples[[t]]
    phis_t = fit$phis[[t]]
    
    c_t = as.integer(c_t)
    mu  = phis_t[c_t]  # mean for each observation under current allocation
    
    ll[t] = sum(dnorm(Y, mean = mu, sd = sd_y, log = TRUE))
  }
  
  ll
}

# ---- helper: relabel cluster ids to 1..H and return mapping ----
# used by dp_slice_mixmodel_normal_normal, 
# dp_CRPwithAtoms_mixmodel_normal_normal, and 
# dp_CRPnoAtoms_mixmodel_normal_normal

relabel_to_consecutive <- function(c) {
  labs   = sort(unique(c))
  new_id = match(c, labs) 
  list(c = new_id, labs = labs)
}