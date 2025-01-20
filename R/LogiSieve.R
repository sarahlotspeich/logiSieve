#' Sieve maximum likelihood estimator (SMLE) for two-phase logistic regression problems with covariate measurement error
#'
#' This function returns the sieve maximum likelihood estimators (SMLE) for the logistic regression model from Lotspeich et al. (2025+).
#'
#' @param analysis_formula formula, analysis model formula (or coercible to formula), a formula expression as for other regression models. The response should be the logistic regression model outcome.
#' @param error_formula formula, covariate error model formula (or coercible to formula), a formula expression as for other regression models. The response should be the error-free version of the error-prone of the covariate, and the covariate should be the names of the B-spline columns.
#' @param data dataframe, a dataframe with one row per subject containing all variables from \code{analysis_formula} and \code{error_formula}.
#' @param initial_lr_params character, initial values for parametric model parameters. Choices include (1) \code{"Zero"} (non-informative starting values) or (2) \code{"Complete-data"} (estimated based on validated subjects only)
#' @param pert_scale scalar, size of the perturbation used in estimating the standard errors via profile likelihood. If none is supplied, default is \code{pert_scale = 1}.
#' @param no_se logical, indicator for whether standard errors are desired. Defaults to \code{no_se = FALSE}.
#' @param tol scalar, tolerance between iterations in the EM algorithm used to define convergence.
#' @param max_iter scalar, maximum number of iterations allowed in the EM algorithm.
#' @param output character, level of fitted model output to be returned. Defaults to \code{output = "logORs"}, but \code{output = "all"} is also possible.
#' @return A list with the following named slots:
#' \item{model_coeff}{dataframe with final model coefficients and standard error estimates (where applicable) for the analysis model.}
#' \item{bspline_coeff}{dataframe with B-spline coefficients for the covariate error model. (Only returned if \code{output = "all"}.)}
#' \item{vcov}{covariance matrix of \code{model_coeff} for the analysis model.}
#' \item{converged}{indicator of EM algorithm convergence for parameter estimates.}
#' \item{se_converged}{indicator of standard error estimate convergence.}
#' \item{converged_msg}{(where applicable) description of non-convergence.}
#' \item{iterations}{number of iterations completed by EM algorithm to find parameter estimates. (Only returned if \code{output = "all"}.)}
#' \item{od_loglik_at_conv}{value of the observed-data log-likelihood at convergence. (Only returned if \code{output = "all"}.)}
#' @export
#' @importFrom stats as.formula
#' @importFrom stats glm

logiSieve = function(analysis_formula, error_formula, data, initial_lr_params = "Zero", 
                     pert_scale = 1, no_se = FALSE, tol = 1E-4, max_iter = 1000, output = "logORs")
{
  # Extract variable names from user-specified formulas
  Y = as.character(as.formula(analysis_formula))[2] ## outcome
  X_val = as.character(as.formula(error_formula))[2] ## error-free covariate
  C = setdiff(x = unlist(strsplit(x = gsub(pattern = " ",
                                           replacement = "",
                                           x = as.character(as.formula(analysis_formula))[3]),
                                  split = "+",
                                  fixed = TRUE)),
              y = X_val)
  Bspline = unlist(strsplit(x = gsub(pattern = " ",
                                     replacement = "",
                                     x = as.character(as.formula(error_formula))[3]),
                            split = "+",
                            fixed = TRUE))
  
  # Create validation indicator 
  data$V = as.numeric(!is.na(data[, X_val])) ## = 1 if validated, = 0 otherwise
  
  # Save sample sizes ---------------------------------------------
  N = nrow(data) ## total (phase I)
  n = sum(data[, "V"]) ## validated (phase II)

  # Reorder so that the n validated subjects are first ------------
  data = data[order(as.numeric(data[, "V"]), decreasing = TRUE), ]

  # Add the B spline basis -------------------------------------------
  sn = ncol(data[, Bspline])
  if(0 %in% colSums(data[c(1:n), Bspline])) {
    warning("Empty sieve in validated data. Reconstruct B-spline basis and try again.", call. = FALSE)
    if(output == "logORs") {
      return(list(model_coeff = data.frame(coeff = NA, se = NA),
                  vcov = NA,
                  converged = NA,
                  se_converged = NA,
                  converged_msg = "B-spline error"))
    } else {
      return(list(model_coeff = data.frame(coeff = NA, se = NA),
                  bspline_coeff = NA,
                  vcov = NA,
                  converged = NA,
                  se_converged = NA,
                  converged_msg = "B-spline error",
                  iterations = 0,
                  od_loglik_at_conv = NA))
    }
  }
  # ------------------------------------------- Add the B spline basis

  # Save distinct X -------------------------------------------------
  x_obs = data.frame(unique(data[1:n, c(X_val)]))
  x_obs = data.frame(x_obs[order(x_obs[, 1]), ])
  m = nrow(x_obs)
  x_obs_stacked = do.call(what = rbind, 
                          args = replicate(n = (N - n), 
                                           expr = x_obs, 
                                           simplify = FALSE))
  x_obs_stacked = data.frame(x_obs_stacked[order(x_obs_stacked[, 1]), ])
  colnames(x_obs) = colnames(x_obs_stacked) = c(X_val)

  # Save static (X*,X,Y,C) since they don't change ---------------
  comp_dat_val = data[c(1:n), c(Y, X_val, C, Bspline)]
  comp_dat_val = merge(x = comp_dat_val, 
                       y = data.frame(x_obs, k = 1:m), 
                       all.x = TRUE)
  comp_dat_val = comp_dat_val[, c(Y, X_val, C, Bspline, "k")]
  comp_dat_val = data.matrix(comp_dat_val)

  # (m x n)xd vectors of each (one column per person, one row per x) --
  comp_dat_unval = suppressWarnings(
    data.matrix(
      cbind(data[-c(1:n), c(Y, C, Bspline)],
            x_obs_stacked,
            k = rep(seq(1, m), each = (N - n)))
    )
  )
  comp_dat_unval = comp_dat_unval[, c(Y, X_val, C, Bspline, "k")]
  comp_dat_all = rbind(comp_dat_val, comp_dat_unval)

  # Initialize B-spline coefficients {p_kj}  ------------
  ## Numerators sum B(Xi*) over k = 1,...,m -------------
  ## Save as p_val_num for updates ----------------------
  ## (contributions don't change) -----------------------
  p_val_num = rowsum(x = comp_dat_val[, Bspline], 
                     group = comp_dat_val[, "k"], 
                     reorder = TRUE)
  prev_p = p0 =  t(t(p_val_num) / colSums(p_val_num))
  
  theta_design_mat = cbind(int = 1, comp_dat_all[, c(X_val, C)])

  # Initialize parameter values -------------------------------------
  if(!(initial_lr_params %in% c("Zero", "Complete-data"))) {
    message("Invalid starting values provided. Non-informative zeros assumed.")
    initial_lr_params = "Zero"
  }
  if(initial_lr_params == "Zero") {
    prev_theta = theta0 = matrix(data = 0, 
                                 nrow = ncol(theta_design_mat), 
                                 ncol = 1)
  } else if(initial_lr_params == "Complete-data") {
    prev_theta = theta0 = matrix(glm(formula = analysis_formula, 
                                     family = "binomial", 
                                     data = data.frame(data[c(1:n), ]))$coefficients, 
                                 ncol = 1)
  }

  # Estimate theta and p using EM -----------------------------------
  CONVERGED = FALSE
  CONVERGED_MSG = "Unknown"
  it = 1
  while(it <= max_iter & !CONVERGED) {
    # E Step ----------------------------------------------------------
    ## Update the psi_kyji for unvalidated subjects -------------------
    ### P(Y|X) --------------------------------------------------------
    mu_theta = as.numeric(theta_design_mat[-c(1:n), ] %*% prev_theta)
    pY_X = 1 / (1 + exp(- mu_theta))
    I_y0 = comp_dat_unval[, Y] == 0
    pY_X[I_y0] = 1 - pY_X[I_y0]
    ### -------------------------------------------------------- P(Y|X)
    ###################################################################
    ### P(X|X*) -------------------------------------------------------
    ### p_kj ----------------------------------------------------------
    ### need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
    ### multiply by the B-spline terms
    pX = prev_p[rep(seq(1, m), each = (N - n)), ] * 
      comp_dat_unval[, Bspline]
    ### ---------------------------------------------------------- p_kj
    ### ------------------------------------------------------- P(X|X*)
    ###################################################################
    ### Estimate conditional expectations -----------------------------
    ### P(Y|X,C)p_kjB(X*) -------------------------------------------
    psi_num = c(pY_X) * pX
    ### Update denominator ------------------------------------------
    #### Sum up all rows per id (e.g. sum over xk) ------------------
    psi_denom = rowsum(psi_num, group = rep(seq(1, (N - n)), times = m))
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
    ## Update theta using weighted logistic regression ----------------
    ### Gradient ------------------------------------------------------
    mu = theta_design_mat %*% prev_theta
    w_t = c(rep(1, n), w_t)
    gradient_theta = matrix(data = c(colSums(w_t * c((comp_dat_all[, Y] - 1 + exp(-mu) / (1 + exp(- mu)))) * theta_design_mat)), ncol = 1)
    ### ------------------------------------------------------ Gradient
    ### Hessian -------------------------------------------------------
    post_multiply = c((exp(- mu) / (1 + exp(- mu))) * (exp(- mu)/(1 + exp(- mu)) - 1)) * w_t * theta_design_mat
    hessian_theta = apply(theta_design_mat, MARGIN = 2, FUN = hessian_row, pm = post_multiply)
    new_theta = tryCatch(expr = prev_theta - solve(hessian_theta) %*% gradient_theta,
                          error = function(err) {
                            matrix(NA, nrow = nrow(prev_theta))
                          })
    if (any(is.na(new_theta))) {
      new_theta = suppressWarnings(matrix(glm(formula = analysis_formula, 
                                              family = "binomial", 
                                              data = data.frame(comp_dat_all), 
                                              weights = w_t)$coefficients, 
                                          ncol = 1))
    }
    ### Check for convergence -----------------------------------------
    theta_conv = abs(new_theta - prev_theta) < tol
    ## --------------------------------------------------- Update theta
    ###################################################################
    ## Update {p_kj} --------------------------------------------------
    ### Update numerators by summing u_t over i = 1, ..., N -----------
    new_p_num = p_val_num +
      rowsum(psi_t, group = rep(seq(1, m), each = (N - n)), reorder = TRUE)
    new_p = t(t(new_p_num) / colSums(new_p_num))
    ### Check for convergence ---------------------------------------
    p_conv = abs(new_p - prev_p) < tol
    ## -------------------------------------------------- Update {p_kj}
    ###################################################################
    # M Step ----------------------------------------------------------
    ###################################################################

    all_conv = c(theta_conv, p_conv)
    if (mean(all_conv) == 1) { CONVERGED = TRUE }

    # Update values for next iteration  -------------------------------
    it = it + 1
    prev_theta = new_theta
    prev_p = new_p
    #  ------------------------------- Update values for next iteration
  }

  rownames(new_theta) = c("Intercept", X_val, C)

  if(!CONVERGED) {
    if(it > max_iter) { CONVERGED_MSG = "max_iter reached" }
    if(output == "logORs") {
      return(list(model_coeff = data.frame(coeff = rep(NA, times = nrow(prev_theta)), se = NA),
                  vcov = matrix(NA, nrow = length(new_theta), ncol = length(new_theta)),
                  converged = CONVERGED,
                  se_converged = NA,
                  converged_msg = CONVERGED_MSG))
    } else {
      return(list(model_coeff = data.frame(coeff = rep(NA, times = nrow(prev_theta)), se = NA),
                  bspline_coeff = cbind(x_obs, NA),
                  vcov = matrix(NA, nrow = length(new_theta), ncol = length(new_theta)),
                  converged = CONVERGED,
                  se_converged = NA,
                  converged_msg = CONVERGED_MSG,
                  iterations = it,
                  od_loglik_at_conv = NA))
    }
  }
  if(CONVERGED) { CONVERGED_MSG = "Converged" }
  # ---------------------------------------------- Estimate theta using EM
  if(no_se) {
    ## Calculate pl(theta) -------------------------------------------------
    od_loglik_theta = observed_data_loglik(N = N,
                                           n = n,
                                           Y = Y,
                                           X_val = X_val,
                                           C = C,
                                           Bspline = Bspline,
                                           comp_dat_all = comp_dat_all,
                                           theta = new_theta,
                                           p = new_p)
    
    if(output == "logORs") { 
      return(list(model_coeff = data.frame(coeff = new_theta,
                                           se = NA),
                  vcov = matrix(data = NA, 
                                nrow = length(new_theta), 
                                ncol = length(new_theta)),
                  converged = CONVERGED,
                  se_converged = NA,
                  converged_msg = CONVERGED_MSG))
    } else {
      return(list(model_coeff = data.frame(coeff = new_theta,
                                           se = NA),
                  bspline_coeff = cbind(x_obs, new_p),
                  vcov = matrix(data = NA, 
                                nrow = length(new_theta), 
                                ncol = length(new_theta)),
                  converged = CONVERGED,
                  se_converged = NA,
                  converged_msg = CONVERGED_MSG,
                  iterations = it,
                  od_loglik_at_conv = od_loglik_theta))
    }
  } else {
    # Estimate Cov(theta) using profile likelihood -------------------------
    h_N = pert_scale * N ^ ( - 1 / 2) # perturbation ----------------------------

    ## Calculate pl(theta) -------------------------------------------------
    od_loglik_theta = observed_data_loglik(N = N,
                                           n = n,
                                           Y = Y,
                                           X_val = X_val,
                                           C = C,
                                           Bspline = Bspline,
                                           comp_dat_all = comp_dat_all,
                                           theta = new_theta,
                                           p = new_p)

    I_theta = matrix(data = od_loglik_theta, 
                     nrow = nrow(new_theta), 
                     ncol = nrow(new_theta))

    single_pert_theta = sapply(X = seq(1, ncol(I_theta)),
                               FUN = pl_theta,
                               theta = new_theta,
                               h_N = h_N,
                               n = n,
                               N = N,
                               pYgivX_unval = pY_X,
                               Y = Y,
                               X_val = X_val,
                               C = C,
                               Bspline = Bspline,
                               comp_dat_all = comp_dat_all,
                               p0 = new_p,
                               p_val_num = p_val_num,
                               tol = tol,
                               max_iter = max_iter)

    if (any(is.na(single_pert_theta))) {
      I_theta = matrix(data = NA, 
                       nrow = nrow(new_theta), 
                       ncol = nrow(new_theta))
      SE_CONVERGED = FALSE
    } else {
      spt_wide = matrix(data = rep(c(single_pert_theta), 
                            times = ncol(I_theta)),
                        ncol = ncol(I_theta),
                        byrow = FALSE)
      #for the each kth row of single_pert_theta add to the kth row / kth column of I_theta
      I_theta = I_theta - spt_wide - t(spt_wide)
      SE_CONVERGED = TRUE
    }

    for (c in 1:ncol(I_theta)) {
      pert_theta = new_theta
      pert_theta[c] = pert_theta[c] + h_N
      double_pert_theta = sapply(X = seq(c, ncol(I_theta)),
                                 FUN = pl_theta,
                                 theta = pert_theta,
                                 h_N = h_N,
                                 n = n,
                                 N = N,
                                 pYgivX_unval = pY_X,
                                 Y = Y,
                                 X_val = X_val,
                                 C = C,
                                 Bspline = Bspline,
                                 comp_dat_all = comp_dat_all,
                                 p0 = new_p,
                                 p_val_num = p_val_num,
                                 max_iter = max_iter,
                                 tol = tol)
      dpt = matrix(data = 0, 
                   nrow = nrow(I_theta), 
                   ncol = ncol(I_theta))
      dpt[c,c] = double_pert_theta[1] #Put double on the diagonal
      if(c < ncol(I_theta)) {
        ## And fill the others in on the cth row/ column
        dpt[c, -(1:c)] = dpt[-(1:c), c] = double_pert_theta[-1]
      }
      I_theta = I_theta + dpt
    }

    I_theta = h_N ^ (- 2) * I_theta

    cov_theta = tryCatch(expr = - solve(I_theta),
                         error = function(err) {
                            matrix(data = NA, 
                                   nrow = nrow(I_theta), 
                                   ncol = ncol(I_theta))
                           }
    )
    # ------------------------- Estimate Cov(theta) using profile likelihood
    # if(any(diag(cov_theta) < 0)) {
    #   warning("Negative variance estimate. Increase the pert_scale parameter and repeat variance estimation.")
    #   SE_CONVERGED = FALSE
    # }
    se_theta = tryCatch(expr = sqrt(diag(cov_theta)),
                            warning = function(w) {
                              matrix(NA, nrow = nrow(prev_theta))
                            }
    )
    if (any(is.na(se_theta))) { SE_CONVERGED = FALSE} else { TRUE }
    if(output == "logORs") { 
      return(list(model_coeff = data.frame(coeff = new_theta,
                                           se = se_theta),
                  vcov = cov_theta,
                  converged = CONVERGED,
                  se_converged = SE_CONVERGED,
                  converged_msg = CONVERGED_MSG))
    } else {
      return(list(model_coeff = data.frame(coeff = new_theta,
                                           se = se_theta),
                  bspline_coeff = cbind(x_obs, new_p),
                  vcov = cov_theta,
                  converged = CONVERGED,
                  se_converged = SE_CONVERGED,
                  converged_msg = CONVERGED_MSG,
                  iterations = it,
                  od_loglik_at_conv = od_loglik_theta))
    }
  }
}