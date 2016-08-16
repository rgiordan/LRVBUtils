NewtonsMethod <- function(EvalFun, EvalGrad, EvalHess, theta_init,
                          max_iters=50, tol=1e-12, fn_scale=1, verbose=FALSE) {
  iter <- 1
  done <- FALSE
  f0 <- fn_scale * EvalFun(theta_init)
  theta <- theta_init
  while(iter <= max_iters && !done) {
    hess <- fn_scale * EvalHess(theta)
    hess_ev <- eigen(hess)
    pos_def <- FALSE
    step_max <- Inf
    initial_step <- 1
    grad <- fn_scale * EvalGrad(theta)
    if (min(hess_ev$values) < -1e-8) {
      # Minimizes a function with line search starting at x in direction v
      if (verbose) cat("\n\n\nSearching along negative eigenvector.  ev = ",
                       min(hess_ev$values), "\n\n")
      min_ind <- which.min(hess_ev$values)
      ev <- hess_ev$vectors[, min_ind]
      step_direction <- sign(sum(-1 * grad * ev)) * ev
      initial_step <- 1
    } else {
      pos_def <- TRUE
      step_direction <- as.numeric(-1 * solve(hess, grad))
      step_max <- 1
    }
    ls_result <- LineSearch(EvalFun, EvalGrad, theta, step_direction,
                            step_scale=0.5, max_iters=5000,
                            step_max=step_max, initial_step=initial_step,
                            fn_scale=fn_scale, verbose=FALSE)
    f1 <- fn_scale * EvalFun(ls_result$x)
    diff <- f1 - f0
    theta <- ls_result$x
    f0 <- f1
    if (!is.numeric(diff)) {
      warn("Non-numeric function evaluation.")
      iter <- max_iters
    }
    if (verbose) cat(" iter: ", iter,
                     " diff: ", diff,
                     " f: ", f1,
                     " step: ", ls_result$step_size,
                     "\n")
    if (pos_def && abs(diff) < tol) {
      if (verbose) cat("Done.\n")
      done <- TRUE
    }
    iter <- iter + 1
  }

  return(list(theta=theta, done=done))
}
