#' Observed-data log-likelihood
#'
#' This function returns the value of the observed-data log-likelihood (equation (2) in Lotspeich et al. (2021))
#' for a given dataset and parameter values `theta` and `p`.
#
#'
#' @param N Phase I sample size
#' @param n Phase II sample size
#' @param Y Column names with the validated outcome.
#' @param X_val Column name(s) with the validated predictors. 
#' @param C (Optional) Column name(s) with additional error-free covariates.
#' @param Bspline Vector of column names containing the B-spline basis functions.
#' @param comp_dat_all Augmented dataset containing rows for each combination of unvalidated subjects' data with values from Phase II (a matrix)
#' @param theta Parameters for the analysis model (a column vector)
#' @param p B-spline coefficients for the approximated covariate error model (a matrix)
#' @return Scalar value of the function
#' @export

observed_data_loglik <- function(N, n, Y = NULL, X_val = NULL, C = NULL, Bspline = NULL, 
                                 comp_dat_all, theta, p) {
  #sn <- ncol(p)
  m <- nrow(p)
  
  # For validated subjects --------------------------------------------------------
  #################################################################################
  ## Sum over log[P_theta(Yi|Xi)] -------------------------------------------------
  pY_X <- 1 / (1 + exp(-as.numeric(cbind(int = 1, comp_dat_all[c(1:n), c(X_val, C)]) %*% theta)))
  pY_X <- ifelse(as.vector(comp_dat_all[c(1:n), c(Y)]) == 0, 1 - pY_X, pY_X)
  return_loglik <- sum(log(pY_X))
  ## ------------------------------------------------- Sum over log[P_theta(Yi|Xi)]
  #################################################################################
  ## Sum over I(Xi=xk)Bj(Xi*)log p_kj ---------------------------------------------
  pX <- p[comp_dat_all[c(1:n), "k"], ]
  log_pX <- log(pX)
  log_pX[log_pX == -Inf] <- 0
  return_loglik <- return_loglik + sum(comp_dat_all[c(1:n), Bspline] * log_pX, na.rm = TRUE)
  ## --------------------------------------------- Sum over I(Xi=xk)Bj(Xi*)log q_kj
  #################################################################################
  # -------------------------------------------------------- For validated subjects

  # For unvalidated subjects ------------------------------------------------------
  ## Calculate P_theta(y|x) for all (y,xk) ----------------------------------------
  pY_X <- 1 / (1 + exp(-as.numeric(cbind(int = 1, comp_dat_all[-c(1:n), c(X_val, C)]) %*% theta)))
  pY_X[which(comp_dat_all[-c(1:n), Y] == 0)] <- 1 - pY_X[which(comp_dat_all[-c(1:n), Y] == 0)]
  ## ---------------------------------------- Calculate P_theta(y|x) for all (y,xk)
  ################################################################################
  ## Calculate Bj(Xi*) p_kj for all (k,j) ----------------------------------------
  pX <- rowSums(p[comp_dat_all[-c(1:n), "k"], ] * comp_dat_all[-c(1:n), Bspline])
  ## ---------------------------------------- Calculate Bj(Xi*) p_kj for all (k,j)
  ################################################################################
  ## Calculate sum of P(y|xk) x Bj(X*) x p_kj ------------------------------------
  person_sum <- rowsum(pY_X * pX, 
                       group = rep(seq(1, (N - n)), 
                                   times = m))
  log_person_sum <- log(person_sum)
  log_person_sum[log_person_sum == -Inf] <- 0
  ## And sum over them all -------------------------------------------------------
  return_loglik <- return_loglik + sum(log_person_sum)
  ################################################################################
  # ----------------------------------------------------- For unvalidated subjects
  return(return_loglik)
}
