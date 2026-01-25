# Slice-Sampler 

# -------------------------------------------------------------------------
# Main sampler
# -------------------------------------------------------------------------
dp_slice_mixmodel_normal_normal <- function(
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
    
  lik  <- function(y, phi){
    dnorm(y, phi, sqrt(sigma2), log = TRUE)
  }
    
  post <- function(y_k, n_k){
    # Posterior of phi_k | y_idx ~ N(mean_n, tau2_n)
    tau2_n  =  1 / (n_k / sigma2 + 1 / tau20)
    mu_n    =  tau2_n * (sum(y_k) / sigma2 + mu0 / tau20)
    rnorm(1, mu_n, sqrt(tau2_n))
  }
    
  prior <- function(x){rnorm(1, mu0, sqrt(tau20))}
  
  if (!is.null(seed)) set.seed(seed)
  
  n = length(Y)
  if (is.null(c_init)) {
    # simple init: all in one cluster
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
  
  pb <- progress_bar$new(
    format = " MCMC [:bar] :percent Estimated completion time: :eta",
    total = Tot, clear = FALSE, width= 100)
  
  for (t in seq_len(Tot)) {
    
    # counts
    n_h = tabulate(c, nbins = H)  # n_h[h] = sum_i 1(c_i=h)
    
    # ---- (2a) sample (pi_1,...,pi_H, pi_star) ~ Dirichlet(n_1,...,n_H, alpha) ----
    alpha_vec = c(n_h, alpha)
    dir_draw  = rdirichlet(n = 1, alpha_vec)
    pi        = dir_draw[seq_len(H)]
    pi_star   = dir_draw[H + 1L]
    # ---- (2b) sample phis from full-conditionals ----
    phis = NULL
    for(h in seq_len(H)){
      phis[h] = post(Y[c==h], n_h[h])
    }
    
    
    # ---- (3) sample slices u_i ~ Unif(0, pi_{c_i}) ----
    u      = runif(n, min = 0, max = pi[c])
    u_star = min(u)
    
    # ---- (4) expand components until pi_star <= u_star ----
    K = H
    while (pi_star > u_star) {
      K = K + 1L
      # sample new atom
      phis[K] = prior(0)
      
      # sample stick piece from remaining mass
      vK      = rbeta(1L, shape1 = 1, shape2 = alpha)
      piK     = vK * pi_star
      pi      = c(pi, piK)
      pi_star = pi_star - piK
    }
    
    # ---- (5) update allocations ----
    for (i in seq_len(n)) {
      Ai = which(pi > u[i]) # active components for obs i
      if(length(Ai) == 1){
        c[i] = Ai
      }else{
        # evaluate log-likelihoods for k in Ai
        logw = lik(Y[[i]], phis[Ai]) 
        m    = max(logw)
        w    = exp(logw - m)
        c[i] = sample(Ai, size = 1, prob = w)
      }
    }
    # ---- (1) H and relabel c to 1..H ----
    rel = relabel_to_consecutive(c)
    c   = rel$c
    H   = length(unique(c))
    
    alpha = rgamma(1, 3 + H, 3*log(n) - log(rbeta(1, 1, alpha, n)) )
    
    # store
    store_idx = store_idx + 1L
    out_c[[store_idx]]    = c
    out_K[store_idx]      = K
    out_H[store_idx]      = H
    out_phis[[store_idx]] = phis
    
    if(t == 10){
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


# -------------------------------------------------------------------------
# Example usage (toy):
# -------------------------------------------------------------------------

# Simulate data
# set.seed(1)
# Y <- rnorm(2000, mean = c(rep(-3,1000), rep(3,1000)), sd = 1)

# Run sampler
# fit_SS = dp_slice_mixmodel_normal_normal(Y, seed = 0)

# Compute log lik at each iteration
# loglik_ss = dp_loglik_trace_mixmodel_normal_normal(Y, fit_SS)
# plot(loglik_ss, type = "l", main = "SS")

#num of clusters 
# H = unlist(lapply(fit_SS$c_samples, function(x) length(unique(x))))
# hist(H)
