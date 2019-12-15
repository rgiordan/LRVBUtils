# ######################################
# # Influence function helpers


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
  log_prior <- GetLogPrior(param_draws)
  mcmc_importance_ratio <- exp(log_prior - log(dens_at_draws$y))

  GetConditionalMeanDiff <- function(g_draws) {
      return(GetEYGivenX(x_draws=param_draws, y_draws=g_draws) - mean(g_draws))
  }

  GetMCMCInfluence <- function(g_draws) {
    conditional_mean_diff <-  GetConditionalMeanDiff(g_draws)
    return(conditional_mean_diff *  mcmc_influence_ratio)
  }


  # Get information about the worst-case perturbation of <param> prior on a
  # component with draws <g_draws>.  We are computing
  #
  # dE[g] / du
  #
  # where the prior on <param> is perturbed by
  # p(param) <- (p(param) + u(param)) / normalizer
  # and where ||u||_2 \le 1
  #
  # Args:
  #  - g_draws: Draws from the variable being influenced.
  #
  # Returns:
  #  - worst_case: |dE[g] / du| for the worst-case u
  #  - worst_u: The u that achieves <worst_case>
  #  - mcmc_influence: The influence function of p(param) on E[g].
  #  - l2_var: The approximate variance of <worst_case>.
  GetMCMCWorstCaseResults <- function(g_draws) {
    mcmc_influence <- GetMCMCInfluence(g_draws)
    mcmc_influence_pos <- mcmc_influence * (mcmc_influence > 0)
    mcmc_influence_neg <- -1 * mcmc_influence * (mcmc_influence < 0)

    # Take the expectation of the influence function squared wrt the prior
    # using importance sampling.
    norm_draws_pos <- (mcmc_influence_pos^2) * mcmc_importance_ratio
    norm_draws_neg <- (mcmc_influence_neg^2) * mcmc_importance_ratio

    l2_pos <- sqrt(mean(norm_draws_pos))
    l2_neg <- sqrt(mean(norm_draws_neg))

    # Get the sample variance with the delta method.
    l2_pos_var <-
        var(norm_draws_pos) / (length(norm_draws_pos) * 4 * mean(norm_draws_pos))
    l2_neg_var <-
        var(norm_draws_neg) / (length(norm_draws_neg) * 4 * mean(norm_draws_pos))

    pos_worse <- (l2_pos > l2_neg)
    worst_case <- max(l2_pos, l2_neg)
    if (pos_worse) {
        l2_var <- l2_pos_var
    } else {
        l2_var <- l2_neg_var
    }
    pos_sign <- 2 * pos_worse - 1
    worst_u <-
        exp(log_prior) * abs(mcmc_influence) *
        (pos_sign * mcmc_influence >= 0) / worst_case

    return(list(worst_case=worst_case,
                worst_u=worst_u,
                mcmc_influence=mcmc_influence,
                l2_var=l2_var))
  }

  GetMCMCWorstCase <- function(g_draws) {
    return(GetMCMCWorstCaseResults(g_draws)$worst_case)
  }

  return(list(GetMCMCInfluence=GetMCMCInfluence,
              GetMCMCWorstCase=GetMCMCWorstCase,
              GetMCMCWorstCaseResults=GetMCMCWorstCaseResults,
              GetConditionalMeanDiff=GetConditionalMeanDiff,
              param_draws=param_draws,
              dens_at_draws=dens_at_draws$y,
              log_prior_at_draws=log_prior))
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
  GetLogQGradTerms,
  GetLogQ,
  GetLogPrior) {

  u_draws <- DrawImportanceSamples(num_draws)

  log_prior <- GetLogPrior(u_draws)
  log_q <- GetLogQ(u_draws)
  importance_lp <- GetImportanceLogProb(u_draws)

  importance_lp_ratio <- log_prior - importance_lp
  influence_lp_ratio <- log_q - log_prior
  log_q_grad_terms <- GetLogQGradTerms(u_draws)

  influence_fun  <- t(log_q_grad_terms) * exp(influence_lp_ratio)

  # Take integrals using importance sampling.
  u_influence_mat <- (influence_fun ^ 2) * exp(importance_lp_ratio)
  u_influence_mat_pos <- ((influence_fun > 0) * influence_fun ^ 2) * exp(importance_lp_ratio)
  u_influence_mat_neg <- ((influence_fun < 0) * influence_fun ^ 2) * exp(importance_lp_ratio)

  u_influence_mat_pos_norm <- sqrt(colMeans(u_influence_mat_pos))
  u_influence_mat_neg_norm <- sqrt(colMeans(u_influence_mat_neg))

  # Get the worst-case perturbation and its magnitude.
  pos_worse <- u_influence_mat_pos_norm > u_influence_mat_neg_norm
  worst_case <- ifelse(pos_worse, u_influence_mat_pos_norm, u_influence_mat_neg_norm)
  GetWorstCasePerturbationColumn <- function(ind) {
      # Choose the positive or negative part according to pos_worse
      pos_sign <- 2 * pos_worse[ind] - 1
      exp(log_prior) * abs(influence_fun[, ind]) *
        (pos_sign * influence_fun[, ind] >= 0) / worst_case[ind]
  }
  worst_case_u <- sapply(1:ncol(influence_fun), GetWorstCasePerturbationColumn)

  return(list(
    u_draws=u_draws,
    influence_fun=influence_fun,
    importance_lp_ratio=importance_lp_ratio,
    influence_lp_ratio=influence_lp_ratio,
    log_prior=log_prior,
    log_q=log_q,
    importance_lp=importance_lp,
    log_q_grad_terms=log_q_grad_terms,
    worst_case=worst_case,
    worst_case_u=worst_case_u))
}
