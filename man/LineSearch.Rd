\name{LineSearch}
\alias{LineSearch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LineSearch(EvalFun, EvalGrad, x, v, initial_step = 1, step_scale = 0.8, max_iters = 5000, step_min = 0, step_max = Inf, fn_decrease = 1e-14, grad_decrease = 0.5, fn_scale = 1, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EvalFun}{
%%     ~~Describe \code{EvalFun} here~~
}
  \item{EvalGrad}{
%%     ~~Describe \code{EvalGrad} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{v}{
%%     ~~Describe \code{v} here~~
}
  \item{initial_step}{
%%     ~~Describe \code{initial_step} here~~
}
  \item{step_scale}{
%%     ~~Describe \code{step_scale} here~~
}
  \item{max_iters}{
%%     ~~Describe \code{max_iters} here~~
}
  \item{step_min}{
%%     ~~Describe \code{step_min} here~~
}
  \item{step_max}{
%%     ~~Describe \code{step_max} here~~
}
  \item{fn_decrease}{
%%     ~~Describe \code{fn_decrease} here~~
}
  \item{grad_decrease}{
%%     ~~Describe \code{grad_decrease} here~~
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
function (EvalFun, EvalGrad, x, v, initial_step = 1, step_scale = 0.8, 
    max_iters = 5000, step_min = 0, step_max = Inf, fn_decrease = 1e-14, 
    grad_decrease = 0.5, fn_scale = 1, verbose = FALSE) 
{
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
        f_condition <- (new_f <= f0 + fn_decrease * step_size * 
            slope0)
        grad_condition <- (abs(new_slope) <= grad_decrease * 
            abs(slope0)) || (abs(new_slope) < 1e-07)
        if (verbose) {
            cat(" iter: ", iter, " step_min: ", step_min, " step_max: ", 
                step_max, " step_size: ", step_size, " f diff: ", 
                f0 - new_f, " grad diff: ", abs(slope0) - abs(new_slope), 
                "\n")
        }
        if (f_condition && grad_condition) {
            done <- TRUE
            if (verbose) 
                cat("Done.\n")
        }
        else if (f_condition) {
            step_min <- step_size
            if (is.infinite(step_max)) {
                step_size <- step_size/step_scale
            }
            else {
                step_size <- step_max - (step_max - step_size) * 
                  step_scale
            }
            if (verbose) 
                cat("Increasing step size to ", step_size, "\n")
        }
        else {
            step_max <- step_size
            step_size <- step_min + (step_size - step_min) * 
                step_scale
            if (verbose) 
                cat("Decreasing step size to ", step_size, "\n")
        }
        if (abs(step_max - step_min) < 1e-08) {
            iter <- max_iters
        }
        iter <- iter + 1
    }
    return(list(x = x_new, f = new_f, grad = new_grad, step_size = step_size, 
        done = done))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
