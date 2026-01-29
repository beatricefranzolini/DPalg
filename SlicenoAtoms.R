################################################################################
## this code contains the MCMC function to run the slice sampler with no atoms
#  for DP mixture of univariate Normals
#  Kernel:   y | mu ~ N(mu, sigma2)
#  Base:     mu ~ N(mu0, tau20)
#  DP:       G ~ DP(alpha, N(mu0, tau20)) with alpha random
################################################################################

# -------------------------------------------------------------------------
# Main sampler
# Inputs:   Y      -> data as a n X 1 vector
#           Tot    -> number of iterations 
#           c_init -> initialization for the partition, 
#                     if NULL is st to k-means solution with 5 clusters
#           seed   -> seed 
#           hyper  -> (kernel variance, base measure mean, base measure variance)
# Outputs:   c_samples -> chain of clustering labels
#            phis      -> chain of atoms
#            K         -> chain of dynamic truncation threshold
#            H         -> chain of number of cluster
#            time      -> wall-clock time needed to run the chain
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Example usage (toy):
# -------------------------------------------------------------------------

# Simulate data
# set.seed(1)
# Y <- rnorm(2000, mean = c(rep(-3,1000), rep(3,1000)), sd = 1)

# Run sampler
# fit_SS = dp_slicenoAtoms_mixmodel_normal_normal(Y, seed = 0)
# -------------------------------------------------------------------------

dp_slicenoAtoms_mixmodel_normal_normal <- function(
    Y,
    Tot = 3000,
    c_init = NULL,
    seed   = NULL,
    hyper = c(1,0,1)
) {
  sigma2 = hyper[1]
  mu0    = hyper[2]
  tau20  = hyper[3]
  alpha  = 1
  
  pred_logdens = function(y, n_k, sum_y) {
    tau2_n = 1 / (n_k / sigma2 + 1 / tau20)
    mu_n   = tau2_n * (sum_y / sigma2 + mu0 / tau20)
    dnorm(y, mean = mu_n, sd = sqrt(sigma2 + tau2_n), log = TRUE)
  }
  
  post_sum <- function(sum_k, n_h){
    # Posterior of phi_k | y_idx ~ N(mean_n, tau2_n)
    tau2_n  =  1 / (n_h / sigma2 + 1 / tau20)
    mu_n    =  tau2_n * (sum_k / sigma2 + mu0 / tau20)
    rnorm(H, mu_n, sqrt(tau2_n))
  }
  
  prior <- function(x){rnorm(1, mu0, sqrt(tau20))}
  
  if (!is.null(seed)) set.seed(seed)
  
  n = length(Y)
  if (is.null(c_init)) {
    c = kmeans(Y, 
               centers = min(5,n))$cluster
  } else {
    stopifnot(length(c_init) == n)
    c = as.integer(c_init)
  }
  
  # storage
  out_c    = vector("list", Tot)
  out_phis = vector("list", Tot)
  out_K    = integer(Tot)
  out_H    = integer(Tot)
  
  store_idx  = 0L
  start.time = Sys.time()             #save starting time 
  
  rel = relabel_to_consecutive(c)
  c   = rel$c
  H   = length(unique(c))
  
  build_stats <- function(c, H) {
    n_k   = tabulate(c, nbins = H)
    sum_k = numeric(H)
    for (k in seq_len(H)) sum_k[k] <- sum(Y[c == k])
    list(n_k = n_k, sum_k = sum_k)
  }
  stats = build_stats(c, H)
  sum_k = stats$sum_k
  n_h   = stats$n_k
  
  pb <- progress_bar$new(
    format = " MCMC [:bar] :percent Estimated completion time: :eta",
    total = Tot, clear = FALSE, width= 100)
  
  for (t in seq_len(Tot)) {
    
    alpha_vec = c(n_h, alpha)
    dir_draw  = rdirichlet(n = 1, alpha_vec)
    pi        = dir_draw[seq_len(H)]
    pi_star   = dir_draw[H + 1L]
    
    u      = runif(n, min = 0, max = pi[c])
    u_star = min(u)
    
    K = H
    while (pi_star > u_star) {
      K = K + 1L
      
      vK      = rbeta(1L, shape1 = 1, shape2 = alpha)
      piK     = vK * pi_star
      pi      = c(pi, piK)
      pi_star = pi_star - piK
      
      n_h[K]   = 0
      sum_k[K] = 0
    }
    
    for (i in seq_len(n)) {
      Ai = which(pi > u[i]) # active components for obs i
      if(length(Ai) == 1){
        c[i] = Ai
      }else{
        # remove i from its cluster
        y_i = Y[i]
        n_h[c[i]]   = n_h[c[i]] - 1L
        sum_k[c[i]] = sum_k[c[i]] - y_i
        # evaluate log-likelihoods for k in Ai
        logw = pred_logdens(y_i, n_h[Ai], sum_k[Ai])
        m    = max(logw)
        w    = exp(logw - m)
        c[i] = sample(Ai, size = 1, prob = w)
        n_h[c[i]]   = n_h[c[i]] + 1L
        sum_k[c[i]] = sum_k[c[i]] + y_i
      }
    }
    
    temp_n   = n_h[unique(c)]
    temp_sum = sum_k[unique(c)]
    
    rel = relabel_to_consecutive(c)
    c   = rel$c
    H   = length(unique(c))
    
    n_h[unique(c)]   = temp_n
    sum_k[unique(c)] = temp_sum
    n_h   = n_h[1:H]
    sum_k = sum_k[1:H]
    
    phis = post_sum(sum_k, n_h)
    
    alpha = rgamma(1, 3 + H, 3*log(n) - log(rbeta(1, 1, alpha, n)) )
    
    # store
    store_idx = store_idx + 1L
    out_c[[store_idx]]    = c
    out_K[store_idx]      = K
    out_H[store_idx]      = H
    out_phis[[store_idx]] = phis
    
    if(t <= 10){
      if(difftime(Sys.time(), start.time, units = "secs")> 1){
        message("1 sec threshold reached - aborted")
      return(list(
        c_samples = out_c, 
        phis      = out_phis,
        K         = out_K,
        H         = out_H,
        time      = Sys.time() - start.time
      ) )
      }
    }
    pb$tick()
    
  }
  end.time <- Sys.time()
  
  return(list(
    c_samples = out_c, 
    phis      = out_phis,
    K         = out_K,
    H         = out_H,
    time      = end.time - start.time
  ) )
}