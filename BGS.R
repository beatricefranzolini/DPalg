# -------------------------------------------------------------------------
# Main sampler
# -------------------------------------------------------------------------
dp_BGS_mixmodel_normal_normal <- function(
    Y,
    L = 10,
    Tot = 3000,
    c_init = NULL,
    seed   = NULL,
    hyper = c(1,0,1)
){
  sigma2 = hyper[1]
  mu0    = hyper[2]
  tau20  = hyper[3]
  alpha  = 1
  
  lik  <- function(y, phi){
    dnorm(y, phi, sqrt(sigma2), log = TRUE)
  }
  
  post <- function(y_k){
    if(n_h[h]>0){
      # Posterior of phi_k | y_idx ~ N(mean_n, tau2_n)
      tau2_n  =  1 / (n_h[h] / sigma2 + 1 / tau20)
      mu_n    =  tau2_n * (sum(y_k) / sigma2 + mu0 / tau20)
      rnorm(1, mu_n, sqrt(tau2_n))
    }else{
      rnorm(1, mu0, sqrt(tau20))
    }
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  n = length(Y)
  if (is.null(c_init)) {
    # simple init: all in one cluster
    c = kmeans(Y, 
               centers = min(5, L))$cluster
  } else {
    stopifnot(length(c_init) == n)
    c = as.integer(c_init)
  }
  phis = NULL
  
  # storage
  out_c    = vector("list", Tot)
  out_phis = vector("list", Tot)
  out_H    = integer(Tot)
  
  store_idx  = 0L
  start.time = Sys.time()             #save starting time 
  
  pb <- progress_bar$new(
    format = " MCMC [:bar] :percent Estimated completion time: :eta",
    total = Tot, clear = FALSE, width= 100)
  
  for (t in seq_len(Tot)) {
    # ---- (1) counts
    n_h = tabulate(c, nbins = L)  # n_h[h] = sum_i 1(c_i=h)
    
    # ---- (2a) sample (pi_1,...,pi_H, pi_star) ~ Dirichlet(n_1,...,n_H, alpha) ----
    pi  = rdirichlet(n = 1, n_h)
    
    # ---- (2b) sample phis from full-conditionals ----
    for(h in seq_len(L)){
      phis[h] = post(Y[c==h])
    }
    
    # ---- (5) update allocations ----
    for (i in seq_len(n)) {
      # evaluate log-likelihoods for k in Ai
      logw = log(pi) + lik(Y[[i]], phis) 
      m    = max(logw)
      w    = exp(logw - m)
      c[i] = sample(L, size = 1, prob = w)
    }
    H = length(unique(c))
    alpha = rgamma(1, 3 + H, 
                   3*log(n) - log(rbeta(1, 1, alpha, n)) )
    
    # store
    store_idx = store_idx + 1L
    out_c[[store_idx]]    = c
    out_phis[[store_idx]] = phis
    out_H[store_idx]      = H
    
    if(t == 10){
      if(difftime(Sys.time(), start.time, units = "secs")> 1){
        message(" 1 sec threshold reached - aborted")
        return(list(
          c_samples = out_c, 
          phis      = out_phis,
          H         = out_H,
          time      = Sys.time() - start.time
        ) )
      }
    }
    pb$tick()
    
  }
  end.time = Sys.time()
  
  return(list(
    c_samples = out_c, 
    phis      = out_phis,
    H         = out_H,
    time      = end.time - start.time
  ) )
}