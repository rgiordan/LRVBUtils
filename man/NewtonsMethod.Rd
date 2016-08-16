\name{NewtonsMethod}
\alias{NewtonsMethod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
NewtonsMethod(EvalFun, EvalGrad, EvalHess, theta_init, max_iters = 20, tol = 1e-12, fn_scale = 1, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EvalFun}{
%%     ~~Describe \code{EvalFun} here~~
}
  \item{EvalGrad}{
%%     ~~Describe \code{EvalGrad} here~~
}
  \item{EvalHess}{
%%     ~~Describe \code{EvalHess} here~~
}
  \item{theta_init}{
%%     ~~Describe \code{theta_init} here~~
}
  \item{max_iters}{
%%     ~~Describe \code{max_iters} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
  \item{fn_scale}{
%%     ~~Describe \code{fn_scale} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (EvalFun, EvalGrad, EvalHess, theta_init, max_iters = 20, 
    tol = 1e-12, fn_scale = 1, verbose = FALSE) 
{
    iter <- 1
    done <- FALSE
    f0 <- fn_scale * EvalFun(theta_init)
    theta <- theta_init
    while (iter <= max_iters && !done) {
        hess <- fn_scale * EvalHess(theta)
        hess_ev <- eigen(hess)
        pos_def <- FALSE
        step_max <- Inf
        initial_step <- 1
        grad <- fn_scale * EvalGrad(theta)
        if (min(hess_ev$values) < -1e-08) {
            if (verbose) 
                cat("\n\n\nSearching along negative eigenvector.  ev = ", 
                  min(hess_ev$values), "\n\n")
            min_ind <- which.min(hess_ev$values)
            ev <- hess_ev$vectors[, min_ind]
            step_direction <- sign(sum(-1 * grad * ev)) * ev
            initial_step <- 1
        }
        else {
            pos_def <- TRUE
            step_direction <- as.numeric(-1 * solve(hess, grad))
            step_max <- 1
        }
        ls_result <- LineSearch(EvalFun, EvalGrad, theta, step_direction, 
            step_scale = 0.5, max_iters = 5000, step_max = step_max, 
            initial_step = initial_step, fn_scale = fn_scale, 
            verbose = FALSE)
        f1 <- fn_scale * EvalFun(ls_result$x)
        diff <- f1 - f0
        theta <- ls_result$x
        f0 <- f1
        if (!is.numeric(diff)) {
            warn("Non-numeric function evaluation.")
            iter <- max_iters
        }
        if (verbose) 
            cat(" iter: ", iter, " diff: ", diff, " f: ", f1, 
                " step: ", ls_result$step_size, "\n")
        if (pos_def && abs(diff) < tol) {
            if (verbose) 
                cat("Done.\n")
            done <- TRUE
        }
        iter <- iter + 1
    }
    return(list(theta = theta, done = done))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
