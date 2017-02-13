######################################
# Influence function helpers

# Get everything necessary for the calculation of influence functions from a draw, u.
GetInfluenceFunctionSampleFunction <- function(
  GetLogVariationalDensity,
  GetLogPrior,
  GetULogDensity,
  lrvb_pre_factor) {

  function(u) {
    log_q_derivs <- GetLogVariationalDensity(u)
    log_prior_val <- GetLogPrior(u)
    log_u_density <- GetULogDensity(u)

    # It doesn't seem to matter if I just use the inverse.
    # lrvb_term <- -1 * lrvb_terms$jac %*% solve(lrvb_terms$elbo_hess, mu_q_derivs$grad)
    lrvb_term <- lrvb_pre_factor %*% log_q_derivs$grad
    influence_function <- exp(log_q_derivs$val - log_prior_val) * lrvb_term

    return(list(u=u,
                lrvb_term=lrvb_term,
                log_q_val=log_q_derivs$val,
                log_prior_val=log_prior_val,
                log_u_density=log_u_density,
                influence_function=influence_function))
  }
}


# Use the output of a draw from GetInfluenceFunctionSampleFunction to get sensitivity
# to a particular contaminating function.
GetSensitivitySampleFunction <- function(GetLogContaminatingPrior) {
  function(u_draw) {
    log_contaminating_prior_val <- GetLogContaminatingPrior(u_draw$u)

    # The vector of sensitivities.
    log_weight <- log_contaminating_prior_val - u_draw$log_u_density
    log_influence_factor <- u_draw$log_q_val - u_draw$log_prior_val
    sensitivity_draw <- exp(log_weight + log_influence_factor) * u_draw$lrvb_term

    # The "mean value theorem" sensitivity
    prior_val <- exp(u_draw$log_prior_val)
    contaminating_prior_val <- exp(log_contaminating_prior_val)
    prior_ratio <- exp(u_draw$log_prior_val - log_contaminating_prior_val)

    mv_sensitivity_draw <-
      exp(log_influence_factor + u_draw$log_prior_val + log_contaminating_prior_val - u_draw$log_u_density) *
      (log_contaminating_prior_val - u_draw$log_prior_val) *
      u_draw$lrvb_term / (contaminating_prior_val - prior_val)

    return(list(sensitivity_draw=sensitivity_draw,
                mv_sensitivity_draw=mv_sensitivity_draw,
                log_weight=log_weight,
                log_influence_factor=log_influence_factor))
  }
}


# Unpack the output of GetSensitivitySampleFunction
UnpackSensitivityList <- function(sensitivity_list) {
    sens_vec_list <- lapply(sensitivity_list, function(entry) { as.numeric(entry$sensitivity_draw) } )
    sens_vec_list_squared <- lapply(sensitivity_list, function(entry) { as.numeric(entry$sensitivity_draw) ^ 2 } )
    sens_vec_mean <- Reduce(`+`, sens_vec_list) / n_samples
    sens_vec_mean_square <- Reduce(`+`, sens_vec_list_squared) / n_samples
    sens_vec_sd <- sqrt(sens_vec_mean_square - sens_vec_mean^2) / sqrt(n_samples)

    mv_sens_vec_list <- lapply(sensitivity_list, function(entry) { as.numeric(entry$mv_sensitivity_draw) } )
    mv_sens_vec_list_squared <- lapply(sensitivity_list, function(entry) { as.numeric(entry$mv_sensitivity_draw) ^ 2 } )
    mv_sens_vec_mean <- Reduce(`+`, mv_sens_vec_list) / n_samples
    mv_sens_vec_mean_square <- Reduce(`+`, mv_sens_vec_list_squared) / n_samples
    mv_sens_vec_sd <- sqrt(mv_sens_vec_mean_square - mv_sens_vec_mean^2) / sqrt(n_samples)

    return(list(sens_vec_mean=sens_vec_mean, sens_vec_sd=sens_vec_sd,
              mv_sens_vec_mean=mv_sens_vec_mean, mv_sens_vec_sd=mv_sens_vec_sd))
}



####################
# New versions:


# Get influence functions for a univariate parameter with draws from its
# posterior, and whose log prior is given by the function GetLogPrior.  Note
# that this relies on the assumption that the prior factorizes.
#
# Args:
# - param_draws: Draws from the posterior of the parameter.
# - GetLogPrior: A function that returns the log prior value at a draw.
#
# Returns:
# - GetMCMCInfluence: A function that takes draws from a parameter <g> and
#    returns the influence function of <param> on <g> evaluated at <param_draws>.
# - GetMCMCWorstCase: A function like GetMCMCInfluence but which returns the
#    L2 worst-case sensitivity.
# - dens_at_draws: The imputed density at the draws. 
GetMCMCInfluenceFunctions <- function(param_draws, GetLogPrior) {
  GetEYGivenX <- function(x_draws, y_draws) {
    e_y_given_x <- loess.smooth(x_draws, y_draws)
    return(approx(e_y_given_x$x, e_y_given_x$y, xout=x_draws)$y)
  }
  
  mcmc_dens <- density(param_draws)
  dens_at_draws <- approx(mcmc_dens$x, mcmc_dens$y, xout=param_draws)
  mcmc_influence_ratio <- exp(log(dens_at_draws$y) - GetLogPrior(param_draws))
  mcmc_importance_ratio <- exp(GetLogPrior(param_draws) - log(dens_at_draws$y))
  
  GetMCMCInfluence <- function(g_draws) {
    conditional_mean_diff <-
        GetEYGivenX(x_draws=param_draws, y_draws=g_draws) - mean(g_draws)
    return(conditional_mean_diff *  mcmc_influence_ratio)
  }

  GetMCMCWorstCase <- function(g_draws) {
    mcmc_influence <- GetMCMCInfluence(g_draws)
    mcmc_influence_pos <- mcmc_influence * (mcmc_influence > 0)
    mcmc_influence_neg <- -1 * mcmc_influence * (mcmc_influence < 0)
    mcmc_wc <- max(
      sqrt(mean((mcmc_influence_pos^2) * mcmc_importance_ratio)),
      sqrt(mean((mcmc_influence_neg^2) * mcmc_importance_ratio)))
    return(mcmc_wc)  
  }
  
  return(list(GetMCMCInfluence=GetMCMCInfluence,
              GetMCMCWorstCase=GetMCMCWorstCase,
              dens_at_draws=dens_at_draws))
}


# Get an importance sampling version of a variational influence function.
# Args:
# - num_draws: The number of draws for the importance sampler
# - DrawImportanceSamples: Take a number of draws and return draws from the
#     importance samplig distribution
# - GetImportanceLogProb: Return the log probability of the importance sampling
#     function at a draw
# - GetLogQGradTerm: For a draw return a numeric matrix that is:
#     -1 * moment_jac %*% solve(vb_fit$hessian, grad_log_q)
# - GetLogQ: For a draw return the log of the variational density
# - GetLogPrior: For a draw return the log prior density 
GetVariationalInfluenceResults <- function(
  num_draws,
  DrawImportanceSamples,
  GetImportanceLogProb,
  GetLogQGradTerm,
  GetLogQ,
  GetLogPrior) {
  
  u_draws <- DrawImportanceSamples(num_draws)

  log_prior <- sapply(u_draws, GetLogPrior)
  log_q <- sapply(u_draws, GetLogQ)
  importance_lp <- sapply(u_draws, GetImportanceLogProb)

  GetImportanceLogRatio <- function(u) { GetLogPrior(u) - GetImportanceLogProb(u) }
  GetInfluenceLogRatio <- function(u) { GetLogQ(u) - GetLogPrior(u) }

  importance_lp_ratio <- log_prior - importance_lp
  influence_lp_ratio <- log_q - log_prior
  influence_fun  <-
    do.call(rbind, lapply(u_draws, function(u) { GetLogQGradTerm(u) })) * exp(influence_lp_ratio)
  u_influence_mat <- (influence_fun ^ 2) * exp(importance_lp_ratio)
  u_influence_mat_pos <- ((influence_fun > 0) * influence_fun ^ 2) * exp(importance_lp_ratio)
  u_influence_mat_neg <- ((influence_fun < 0) * influence_fun ^ 2) * exp(importance_lp_ratio)
  
  worst_case <-
    sapply(1:ncol(influence_fun),
           function(ind) { sqrt(max(mean(u_influence_mat_pos[, ind]),
                                    mean(u_influence_mat_neg[, ind]))) })
  
  return(list(
    u_draws=u_draws,
    influence_fun=influence_fun,
    importance_lp_ratio=importance_lp_ratio,
    influence_lp_ratio=influence_lp_ratio,
    log_prior=log_prior,
    log_q=log_q,
    importance_lp=importance_lp))
}
