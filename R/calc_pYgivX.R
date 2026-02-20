calc_pYgivX = function(data, successes, failures = NULL, theta, analysis_link) {
  mu_theta = as.numeric(data %*% theta) ## beta0 + beta1x + ... 
  ## Number of trials 
  if (length(failures) == 0) {
    trials = rep(1, nrow(data)) ## If Bernoulli, define as 1
  } else {
    trials = successes + failures ## Otherwise, sum successes + failures
  }
  ## Calculate probability of success on a single trial
  if (analysis_link == "logit") {
    pi_theta = 1 / (1 + exp(- mu_theta)) 
  } else if (analysis_link == "log") {
    pi_theta = exp(mu_theta) 
  }
  ## Evaluate binomial PMF
  pYgivX = dbinom(x = successes, 
                  size = trials, 
                  prob = pi_theta)
}