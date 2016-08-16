LineSearch <- function(EvalFun, EvalGrad, x, v,
                       initial_step=1, step_scale=0.8, max_iters=5000,
                       step_min=0, step_max=Inf,
                       fn_decrease=1e-14, grad_decrease=0.5,
                       fn_scale=1, verbose=FALSE) {
  step_size <- initial_step

  f0 <- EvalFun(x) * fn_scale
  grad0 <- EvalGrad(x) * fn_scale
  slope0 <- sum(grad0 * v)

  iter <- 1
  done <- FALSE
  while (iter <= max_iters && !done) {
    x_new <- x + step_size * v
    new_f <- EvalFun(x_new) * fn_scale
    new_grad <- EvalGrad(x_new) * fn_scale
    new_slope <- sum(new_grad * v)

    # Strong Wolfe conditions
    f_condition <- (new_f <= f0 + fn_decrease * step_size * slope0)
    grad_condition <-
      (abs(new_slope) <= grad_decrease * abs(slope0)) || (abs(new_slope) < 1e-7)

    if (verbose) {
      cat(" iter: ", iter,
          " step_min: ", step_min,
          " step_max: ", step_max,
          " step_size: ", step_size,
          " f diff: ", f0 - new_f,
          " grad diff: ", abs(slope0) - abs(new_slope),
          "\n")
    }

    if (f_condition && grad_condition) {
      done <- TRUE
      if (verbose) cat("Done.\n")
    } else if (f_condition) {
      # f decrease but not slope decrease.
      # Increase step size to get a slope decrease
      step_min <- step_size
      if (is.infinite(step_max)) {
        step_size <- step_size/ step_scale
      } else {
        step_size <- step_max - (step_max - step_size) * step_scale
      }
      if (verbose) cat("Increasing step size to ", step_size, "\n")
    } else {
      # Decrease step size to get a function decrease
      step_max <- step_size
      step_size <- step_min + (step_size - step_min) * step_scale
      if (verbose) cat("Decreasing step size to ", step_size, "\n")
    }

    if (abs(step_max - step_min) < 1e-8) {
      iter <- max_iters
    }

    iter <- iter + 1
  }

  return(list(x=x_new, f=new_f, grad=new_grad, step_size=step_size, done=done))
}
