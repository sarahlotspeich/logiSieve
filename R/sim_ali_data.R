#' Generate error-prone, partially validated allostatic load index (ALI) and healthcare utilization data
#'
#' @param tpr true positive rate for the error-prone ALI. Default is \code{tpr = 0.9}.
#' @param fpr false positive rate for the error-prone ALI. Default is \code{fpr = 0.1}.
#' @param audit_recovery proportion of missing data recovered through the validation study. Default is \code{audit_recovery = 1}.
#' @param N total sample size (phase I) for the error-prone EHR data. Default is \code{N = 1000}.
#' @param n validation sample size (phase II) for the chart review data (must have \code{N > n}). Default is \code{n = 100}.
#' @param lambda_age mean of the Poisson distribution for age. Default is \code{lambda_age = 4.566}.
#' @param betas vector of 3 coefficients for the analysis model (intercept, log odds ratio for X, and log odds ratio for Z). Default is \code{betas = c(-1.93, 1.88, 0.10)}.
#' @param pS vector of 10 probabilities for the components of the ALI being "yes". Default is \code{pS = c(0.2500000, 0.9870130, 0.4549098, 0.1450000, 0.0580000, 0.2490119, 0.3138501, 0.3316391, 0.3111111, 0.0000000)}. 
#' @param pM vector of 10 probabilities for the components of the ALI being missing. Default is \code{pM = c(0.996, 0.153, 0.002, 0.000, 0.000, 0.494, 0.213, 0.213, 0.955, 0.983)}. 
#' @return dataframe with the following columns:
#' \item{id}{Simulated identifiers}
#' \item{X}{True ALI}
#' \item{Xstar}{Unvalidated ALI (from the EHR)}
#' \item{Xval}{Validated ALI (from the chart review)}
#' \item{Y}{Indicator of healthcare utilization}
#' \item{Z}{Age at first encounter}
#' @export

sim_ali_data = function(tpr = 0.9, fpr = 0.1, audit_recovery = 1, N = 1000, n = 100,
                        lambda_age = 4.566, betas = c(-1.93, 1.88, 0.10), 
                        pS = c(0.2500000, 0.9870130, 0.4549098, 0.1450000, 0.0580000,
                               0.2490119, 0.3138501, 0.3316391, 0.3111111, 0.0000000),
                        pM = c(0.996, 0.153, 0.002, 0.000, 0.000,
                               0.494, 0.213, 0.213, 0.955, 0.983)) {
  ## Separate beta coefficients
  beta0 = betas[1] ### intercept 
  beta1 = betas[2] ## log OR on X (given Z)
  beta2 = betas[3] ## log OR on Z (given X)
  
  ## Simulate continuous error-free covariate: age at first encounter
  ### from Poisson(lambda_age)
  Z = rpois(n = N,
            lambda = lambda_age)

  ## Simulate error-free (validated) version of error-prone covariate: ALI
  ### Begin with stress indicators (50 per person) from Bernoulli (pS)
  S = rbinom(n = N * 10,
             size = 1,
             prob = rep(x = pS, times = N))
  S_mat = matrix(data = S,
                 nrow = N,
                 ncol = 10,
                 byrow = TRUE)
  X = rowMeans(S_mat)

  ## Simulate error-prone (EHR) version of the covariate: ALI*
  ### Begin with stress indicators (10 per person) from Bernoulli (1 / (1 + exp(-(gamma0 + gamma1 S))))
  gamma0 = - log((1 - fpr) / fpr) ### define intercept such that P(S* = 0|S = 0) = FPR
  gamma1 = - log((1 - tpr) / tpr) - gamma0 ### define slope such that P(S* = 1|S = 1) = TPR
  Sstar = rbinom(n = N * 10,
                 size = 1,
                 prob = 1 / (1 + exp(- (gamma0 + gamma1 * S))))

  ### Simulate missingness in EHR version of the covariate
  Sstar_miss = rbinom(n = N * 10,
                      size = 1,
                      prob = rep(x = pM, times = N))
  Sstar[which(Sstar_miss == 1)] = NA
  Sstar_mat = matrix(data = Sstar,
                     nrow = N,
                     ncol = 10,
                     byrow = TRUE)
  Xstar = rowMeans(Sstar_mat,
                   na.rm = TRUE)

  ## Simulate outcome: healthcare utilization
  Y = rbinom(n = N,
             size = 1,
             prob = (1 + exp(-(beta0 + beta1 * X - beta2 * Z))) ^ (- 1))

  ## Simulate imperfect audit recovery
  recovered = sample(x = which(Sstar_miss == 1),
                     size = audit_recovery * length(which(Sstar_miss == 1)),
                     replace = FALSE)
  not_recovered = setdiff(x = which(Sstar_miss == 1),
                          y = recovered)
  S[not_recovered] = NA ## components not recovered in audit
  S_mat = matrix(data = S,
                 nrow = N,
                 ncol = 10,
                 byrow = TRUE)
  Xval = rowMeans(S_mat,
                  na.rm = TRUE) ## ALI based on recovered components

  ## Create dataset
  dat = data.frame(id = 1:N, X, Xstar, Xval, Y, Z)

  # Return dataset
  return(dat)
}
