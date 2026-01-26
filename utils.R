# Compute log lik at each iteration
dp_loglik_trace_mixmodel_normal_normal <- function(Y, fit, hyper = c(1, 0, 1)) {
  sd_y   = sqrt(hyper[1])
  
  # Make sure Y is a plain numeric vector
  if (is.list(Y)) Y <- unlist(Y)
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
# used by slice and CRF
relabel_to_consecutive <- function(c) {
  labs   = sort(unique(c))
  new_id = match(c, labs) 
  list(c = new_id, labs = labs)
}

# --- helper to make a data frame for ggplot ---
mk_df <- function(x) {
  data.frame(time = seq_along(x), value = as.numeric(x))
}

# --- helper to evaluate rand index ---
rand_index_evaluation <- function(c_true, fit) {
  
  S = length(fit$c_samples)
  if (is.null(fit$phis[[S]])){return(NULL)}
  
  ll = numeric(S)
  
  for (t in seq_len(S)) {
    c_t    = fit$c_samples[[t]]

    c_t = as.integer(c_t)
    
    ll[t] = rand.index(c_true,c_t)
  }
  
  ll
}