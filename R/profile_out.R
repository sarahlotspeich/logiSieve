#' Profiles out nuisance parameters from the observed-data log-likelihood for a given value of theta
#'
#' For a given vector `theta` to parameterize P(Y|X,C), this function repeats the EM algorithm to find
#' the values of `p` at convergence. The resulting parameters are used to find the profile
#' log-likelihood for `theta` by plugging them into the observed-data log-likelihood.
#' This function is used by `pl_theta()`.
#
#' @param theta Parameters for the analysis model (a column vector)
#' @param N Phase I sample size
#' @param n Phase II sample size
#' @param theta_design_mat Design matrix for model P(Y|X) for validated and unvalidated rows. 
#' @param Y Column names with the validated outcome.
#' @param Bspline Vector of column names containing the B-spline basis functions.
#' @param comp_dat_all Augmented dataset containing rows for validated subjects and each combination of unvalidated subjects' data with values from Phase II
#' @param p0 Starting values for `p`, the B-spline coefficients for the approximated covariate error model (a matrix)
#' @param p_val_num Contributions of validated subjects to the numerator for `p`, which are fixed (a matrix)
#' @param tol tolerance between iterations in the EM algorithm used to define convergence.
#' @param max_iter Maximum number of iterations allowed in the EM algorithm.
#' @return Profile likelihood for `theta`: the value of the observed-data log-likelihood after profiling out other parameters.

profile_out <- function(theta, N, n, theta_design_mat, Y, Bspline, 
                        comp_dat_all, p0, p_val_num, tol, max_iter) {
  # Save useful constants -------------------------------------------
  ## Dimensions and starting values ---------------------------------
  sn <- ncol(p0)
  m <- nrow(p0)
  prev_p <- p0
  
  ## Calculate P(Y|X) for theta, since it won't update --------------
  mu_theta = as.numeric(theta_design_mat[-c(1:n), ] %*% theta)
  pY_X = 1 / (1 + exp(- mu_theta))
  I_y0 = comp_dat_all[-c(1:n), Y] == 0
  pY_X[I_y0] = 1 - pY_X[I_y0]
  
  # Estimate p using EM -----------------------------------------------
  CONVERGED = FALSE
  CONVERGED_MSG = "Unknown"
  it = 1
  while(it <= max_iter & !CONVERGED) {
    # E Step ----------------------------------------------------------
    ## Update the psi_kyji for unvalidated subjects -------------------
    ### P(X|X*) -------------------------------------------------------
    ### p_kj ----------------------------------------------------------
    ### need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
    ### multiply by the B-spline terms
    pX = prev_p[rep(seq(1, m), each = (N - n)), ] * 
      comp_dat_all[-c(1:n), Bspline]
    ### ---------------------------------------------------------- p_kj
    ### ------------------------------------------------------- P(X|X*)
    ###################################################################
    ### Estimate conditional expectations -----------------------------
    ### P(Y|X,C)p_kjB(X*) -------------------------------------------
    psi_num = c(pY_X) * pX
    ### Update denominator ------------------------------------------
    #### Sum up all rows per id (e.g. sum over xk) ------------------
    psi_denom = rowsum(x = psi_num, 
                       group = rep(seq(1, (N - n)), times = m))
    #### Then sum over the sn splines -------------------------------
    psi_denom = rowSums(psi_denom)
    #### Avoid NaN resulting from dividing by 0 ---------------------
    psi_denom[psi_denom == 0] = 1
    ### And divide them! --------------------------------------------
    psi_t = psi_num / psi_denom
    ### Update the w_kyi for unvalidated subjects -------------------
    ### by summing across the splines/ columns of psi_t -------------
    w_t = rowSums(psi_t)
    ### ----------------------------- Estimate conditional expectations
    # ---------------------------------------------------------- E Step
    ###################################################################
    
    ###################################################################
    # M Step ----------------------------------------------------------
    ###################################################################
    ## Update {p_kj} --------------------------------------------------
    ### Update numerators by summing u_t over i = 1, ..., N -----------
    new_p_num = p_val_num +
      rowsum(x = psi_t, 
             group = rep(seq(1, m), each = (N - n)), 
             reorder = TRUE)
    new_p = t(t(new_p_num) / colSums(new_p_num))
    ### Check for convergence ---------------------------------------
    p_conv = abs(new_p - prev_p) < tol
    ## -------------------------------------------------- Update {p_kj}
    ###################################################################
    # M Step ----------------------------------------------------------
    ###################################################################
    
    if (mean(p_conv) == 1) { CONVERGED = TRUE }
    # Update values for next iteration  -------------------------------
    it = it + 1
    prev_p = new_p
    #  ------------------------------- Update values for next iteration
  }

  if(it > max_iter & !CONVERGED) {
    CONVERGED_MSG <- "max_iter reached"
    new_p <- matrix(data = NA, 
                    nrow = nrow(p0), 
                    ncol = ncol(p0))
  }
  if(CONVERGED) CONVERGED_MSG <- "converged"
  # ---------------------------------------------- Estimate theta using EM
  return(list("psi_at_conv" = psi_t,
              "p_at_conv" = new_p,
              "converged" = CONVERGED,
              "converged_msg" = CONVERGED_MSG))
}
