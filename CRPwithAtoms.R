################################################################################
## this code contains the MCMC function to run the CRP with atoms
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
# fit_SS = dp_CRPwithAtoms_mixmodel_normal_normal(Y, seed = 0)
# -------------------------------------------------------------------------

dp_CRPwithAtoms_mixmodel_normal_normal <- function(
    Y,
    Tot = 3000,
    c_init = NULL,
    seed = NULL,
    hyper = c(sigma2 = 1, mu0 = 0, tau20 = 1)
) {
  if (!is.null(seed)) set.seed(seed)
  
  Y = as.numeric(Y)
  n = length(Y)
  
  sigma2 = hyper[["sigma2"]]
  mu0    = hyper[["mu0"]]
  tau20  = hyper[["tau20"]]
  alpha  = 1

  # --- prior predictive (new cluster, n_k=0) ---
  # y | new: N(mu0, sigma2 + tau20)
  prior_pred_logdens = function(y) {
    dnorm(y, mean = mu0, sd = sqrt(sigma2 + tau20), log = TRUE)
  }
  
  # --- initialize allocations ---
  if (is.null(c_init)) {
    c = kmeans(Y, centers = min(5,n))$cluster
  } else {
    stopifnot(length(c_init) == n)
    c = as.integer(c_init)
  }
  
  tmp = relabel_to_consecutive(c)
  c   = tmp$c
  
  # --- cluster sufficient statistics ---
  #   n_k: counts per cluster
  #   sum_k: sum of Y in cluster
  build_stats <- function(c) {
    H     = max(c)
    n_k   = tabulate(c, nbins = H)
    sum_k = numeric(H)
    for (k in seq_len(H)) sum_k[k] <- sum(Y[c == k])
    list(n_k = n_k, sum_k = sum_k)
  }
  stats = build_stats(c)
  n_k   = stats$n_k
  sum_k = stats$sum_k
  H     = max(c)
  mus   = numeric(H)
  
  for (k in seq_len(H)) {
    tau2_n = 1 / (n_k[k] / sigma2 + 1 / tau20)
    mu_n   = tau2_n * (sum_k[k] / sigma2 + mu0 / tau20)
    mus[k] = rnorm(1, mean = mu_n, sd = sqrt(tau2_n))
  }
  
  # --- storage ---
  c_samples  = vector("list", Tot)
  H_trace    = integer(Tot)
  mu_samples = vector("list", Tot)
  
  store_idx = 0L
  start.time = Sys.time()             #save starting time 
  
  pb <- progress_bar$new(
    format = " MCMC [:bar] :percent Estimated completion time: :eta",
    total = Tot, clear = FALSE, width= 100)
  
  # --- main Gibbs loop ---
  for (t in seq_len(Tot)) {
    
    for (i in seq_len(n)) {
      k_old = c[i]
      y_i   = Y[i]
      
      # remove i from its cluster
      n_k[k_old]   = n_k[k_old] - 1L
      sum_k[k_old] = sum_k[k_old] - y_i
      
      # if cluster becomes empty, delete it (swap with last cluster for O(1))
      if (n_k[k_old] == 0L) {
        H = length(n_k)
        if (k_old != H) {
          # move last cluster into position k_old
          n_k[k_old]   = n_k[H]
          sum_k[k_old] = sum_k[H]
          mus[k_old]   = mus[H]
          # relabel assignments pointing to H
          c[c == H] = k_old
        }
        n_k   = n_k[-H]
        sum_k = sum_k[-H]
        mus   = mus[-H]
      }
      
      H = length(n_k)
      
      # compute log weights for existing clusters + new cluster
      # existing: proportional to n_k * pred(y_i | cluster k)
      logw_exist <- numeric(H)
      for (k in seq_len(H)) {
        logw_exist[k] <- log(n_k[k]) + dnorm(y_i, mus[k], 1, log = TRUE) 
      }
      
      # new: proportional to alpha * prior_pred(y_i)
      logw_new = log(alpha) + prior_pred_logdens(y_i)
      
      logw = c(logw_exist, logw_new)
      
      # sample new assignment from categorical with log-weights
      m = max(logw)
      w = exp(logw - m)
      w = w / sum(w)
      k_new = sample.int(H + 1L, size = 1L, prob = w)
      
      if (k_new <= H) {
        # assign to existing cluster
        c[i]         = k_new
        n_k[k_new]   = n_k[k_new] + 1L
        sum_k[k_new] = sum_k[k_new] + y_i
      } else {
        # create new cluster at position H+1
        c[i]  = H + 1L
        n_k   = c(n_k, 1L)
        sum_k = c(sum_k, y_i)
        tau2_n = 1 / (1 / sigma2 + 1 / tau20)
        mu_n   = tau2_n * (y_i / sigma2 + mu0 / tau20)
        mus    = c(mus, rnorm(1, mean = mu_n, sd = sqrt(tau2_n)) )
      }
    }
    
    for (k in seq_len(H)) {
      tau2_n = 1 / (n_k[k] / sigma2 + 1 / tau20)
      mu_n   = tau2_n * (sum_k[k] / sigma2 + mu0 / tau20)
      mus[k] = rnorm(1, mean = mu_n, sd = sqrt(tau2_n))
    }
    
    # --- store ---
    store_idx = store_idx + 1L
    c_samples[[store_idx]] = c
    H_trace[store_idx] = length(n_k)
    mu_samples[[store_idx]] <- mus
    
    alpha = rgamma(1, 3 + H, 3*log(n) - log(rbeta(1, 1, alpha, n)) )
    
    if(t <= 10){
      if(difftime(Sys.time(), start.time, units = "secs")> 1){
        message("1 sec threshold reached - aborted")
      return(list(
        c_samples = c_samples,
        H    = H_trace,
        phis = mu_samples,
        time      = Sys.time() - start.time
      ) )
      }
    }
    pb$tick()
    
  }
  end.time <- Sys.time() 
  
  list(
    c_samples = c_samples,
    phis = mu_samples,
    H    = H_trace,
    time = end.time - start.time
  )
}