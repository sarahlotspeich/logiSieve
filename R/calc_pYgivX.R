calc_pYgivX = function(data, success_col, fail_col, theta) {
  mu_theta = as.numeric(data %*% theta) ## beta0 + beta1x + ... 
  if (link == "logit") {
    pi_theta = 1 / (1 + exp(- mu_theta)) ## probability of success on a single trial
    trials = rep(1, nrow(data)) ## number of trials (1 = Bernoulli)
  } else if (link == "log") {
    pi_theta = exp(mu_theta) ## probability of success on a single trial
    trials = data[, success_col] + data[, fail_col] ## number of trials 
  }
  pYgivX = dbinom(x = data[, success_col], 
                  size = trials, 
                  prob = pi_theta)
}