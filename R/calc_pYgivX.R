calc_pYgivX = function(data, successes, failures, theta) {
  mu_theta = as.numeric(data %*% theta) ## beta0 + beta1x + ... 
  if (link == "logit") {
    pi_theta = 1 / (1 + exp(- mu_theta)) ## probability of success on a single trial
    trials = rep(1, nrow(data)) ## number of trials (1 = Bernoulli)
  } else if (link == "log") {
    pi_theta = exp(mu_theta) ## probability of success on a single trial
    trials = successes + failures ## number of trials 
  }
  pYgivX = dbinom(x = successes, 
                  size = trials, 
                  prob = pi_theta)
}