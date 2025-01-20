#' Profile likelihood for theta, the analysis model parameters
#'
#' This function returns the value of the profile log-likelihood for parameters `theta` of the
#' analysis model P(Y|X,C) after perturbing element `k` of `theta` by some small amount `h_N`.
#
#' @param k A numeric index between 1 and the dimension of theta for the element of theta to be perturbed
#' @param theta Parameters for the analysis model (a column vector) at convergence, resulting from the EM algorithm
#' @param h_N Size of the small perturbation in `theta[k]`, by default chosen to be `h_N =  N ^ ( - 1 / 2)`
#' @param N Phase I sample size
#' @param n Phase II sample size
#' @param pY_X P(Y|X) for unvalidated rows at convergence for \code{theta}.
#' @param Bspline Vector of column names containing the B-spline basis functions.
#' @param comp_dat_all Augmented dataset containing rows for each combination of unvalidated subjects' data with values from Phase II (a matrix)
#' @param p0 Starting values for `p`, the B-spline coefficients for the approximated covariate error model (a matrix)
#' @param p_val_num Contributions of validated subjects to the numerator for `p`, which are fixed (a matrix)
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#' @return Profile likelihood for `theta` after perturbing element `k` by `h_N`.

pl_theta <- function(k, theta, h_N, n, N, pY_X, Bspline, comp_dat_all, 
                     p0 = NULL, p_val_num = NULL, TOL, MAX_ITER) {
  pert <- theta
  pert[k] <- pert[k] + h_N
  pl_params <- profile_out(theta = pert,
                           n = n,
                           N = N,
                           pY_X = pY_X,
                           Bspline = Bspline,
                           comp_dat_all = comp_dat_all,
                           p0 = p0,
                           p_val_num = p_val_num,
                           TOL = TOL,
                           MAX_ITER = MAX_ITER)
  if(pl_params$converged) {
    od_loglik_pert <- observed_data_loglik(N = N,
                                           n = n,
                                           Y = Y,
                                           X_val = X_val,
                                           C = C,
                                           Bspline = Bspline,
                                           comp_dat_all = comp_dat_all,
                                           theta = pert,
                                           p = pl_params$p_at_conv)

  } else { od_loglik_pert <- NA }
  return(od_loglik_pert)
}
