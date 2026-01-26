# Collapsed Gibbs sampler (CRP marginal) for DP mixture of univariate Normals
# Kernel:   y | mu ~ N(mu, sigma2)
# Base:     mu ~ N(mu0, tau20)
# DP:       G ~ DP(alpha, N(mu0, tau20))
#
# Updates:  resample each c_i given others, integrating out mus
#
# Returns:  c_samples, H_trace, and optionally sampled mus

dp_CRPnoAtoms_mixmodel_normal_normal <- function(
    Y,
    Tot = 3000,
    hyper = c(sigma2 = 1, mu0 = 0, tau20 = 1),
    c_init = NULL,
    seed = NULL,
    sample_mus = TRUE
) {
  if (!is.null(seed)) set.seed(seed)
  
  Y = as.numeric(Y)
  n = length(Y)
  
  sigma2 = hyper[["sigma2"]]
  mu0    = hyper[["mu0"]]
  tau20  = hyper[["tau20"]]
  alpha  = 1
  
  # --- predictive density under Normal-Normal (mu integrated out) ---
  # For cluster with n_k, sum_y, posterior mu|data: N(mu_n, tau2_n)
  # predictive y_new | data: N(mu_n, sigma2 + tau2_n)
  pred_logdens = function(y, n_k, sum_y) {
    tau2_n = 1 / (n_k / sigma2 + 1 / tau20)
    mu_n   = tau2_n * (sum_y / sigma2 + mu0 / tau20)
    dnorm(y, mean = mu_n, sd = sqrt(sigma2 + tau2_n), log = TRUE)
  }
  
  # --- prior predictive (new cluster, n_k=0) ---
  # y | new: N(mu0, sigma2 + tau20)
  prior_pred_logdens = function(y) {
    dnorm(y, mean = mu0, sd = sqrt(sigma2 + tau20), log = TRUE)
  }
  
  # --- initialize allocations ---
  if (is.null(c_init)) {
    # simple init: random labels among a few clusters
    c = kmeans(Y, centers = min(5,n))$cluster
  } else {
    stopifnot(length(c_init) == n)
    c = as.integer(c_init)
  }
  
  tmp = relabel_to_consecutive(c)
  c   = tmp$c
  
  # --- cluster sufficient statistics ---
  # We maintain:
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
  
  # --- storage ---
  c_samples  = vector("list", Tot)
  H_trace    = integer(Tot)
  mu_samples = if (sample_mus) vector("list", Tot) else NULL
  
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
          # relabel assignments pointing to H
          c[c == H] = k_old
        }
        n_k   = n_k[-H]
        sum_k = sum_k[-H]
      }
      
      H = length(n_k)
      
      # compute log weights for existing clusters + new cluster
      # existing: proportional to n_k * pred(y_i | cluster k)
      logw_exist = numeric(H)
      for (k in seq_len(H)) {
        logw_exist[k] = log(n_k[k]) + pred_logdens(y_i, n_k[k], sum_k[k])
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
      }
    }

    # --- store ---
    store_idx = store_idx + 1L
    c_samples[[store_idx]] = c
    H_trace[store_idx] = length(n_k)
    
    if (sample_mus) {
      # sample cluster means from posterior given current clusters (optional)
      H   = length(n_k)
      mus = numeric(H)
      for (k in seq_len(H)) {
        tau2_n = 1 / (n_k[k] / sigma2 + 1 / tau20)
        mu_n   = tau2_n * (sum_k[k] / sigma2 + mu0 / tau20)
        mus[k] = rnorm(1, mean = mu_n, sd = sqrt(tau2_n))
      }
      mu_samples[[store_idx]] <- mus
    }
    alpha = rgamma(1, 3 + H, 3*log(n) - log(rbeta(1, 1, alpha, n)) )
    
    if(t <= 10){
      if(difftime(Sys.time(), start.time, units = "secs")> 1){
        message("1 sec threshold reached - aborted")
      return(list(
        c_samples = c_samples,
        H    = H_trace,
        phis = mu_samples,
        time = Sys.time() - start.time
      ) )
      }
    }  
    pb$tick()
  }
  
  end.time <- Sys.time() 
  list(
    c_samples = c_samples,
    H    = H_trace,
    phis = mu_samples,
    time = end.time - start.time
  )
}

# -------------------------
# Example
# -------------------------
# Simulate data
# set.seed(1)
# Y <- rnorm(1000, mean = c(rep(-3,500), rep(3,500)), sd = 1)

# fit_crp <- dp_CRPnoAtoms_mixmodel_normal_normal(Y, seed = 0)
# summary(fit_crp$H)

# Compute log lik at each iteration
# loglik_crp = dp_loglik_trace_mixmodel_normal_normal(Y, fit_crp)
# plot(loglik_crp, type = "l", main = "CRP")

# hist(fit_crp$H)

