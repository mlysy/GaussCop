#' Generalized Box-Cox transformation.
#'
#' @param x Vector of quantiles at which to compute the transformation.
#' @param lambda Exponent of the transformation.  See Details.
#' @param alpha Offset of the transformation.  See Details.
#' @param normalize Logical; if \code{TRUE} divides by the geometric mean.  See Details.
#' @param jacobian Logical; if \code{TRUE} calculates the Jacobian \code{|dz / dx|}, which converts transformed density values back to the original scale.
#' @details The Generalized Power or Box-Cox transformation is
#' \deqn{
#' z = \begin{array}{rl} ((x + \alpha)^\lambda - 1) / (\lambda C^{\lambda-1}) & \lambda \neq 0 \\ C \log(x + \alpha) & \lambda = 0, \end{array}
#' }{
#' z = \begin{array}{rl} ((x + \alpha)^\lambda - 1) / (\lambda C^{\lambda-1}) & \lambda \neq 0 \\ C \log(x + \alpha) & \lambda = 0, \end{array}
#' }
#'
#' where \eqn{C}{C} is the Geometric mean, i.e., \code{C = exp(mean(log(x + alpha)))}.  Note that \code{C} is only calculated if \code{normalize = TRUE}.
#' @return The vector \code{z} of transformed values, and optionally the Jacobian of the inverse transformation.  See Details.
#' @examples
#' # generate data and plot
#' # apply power transform and superimpose on plot
#' # finally, superimpose N(0, 1) on plot
#' n <- 1e5
#' df <- 5
#' X <- rchisq(n, df = df)
#' xdens <- kernelXD(X)
#' xdens.trans <- kernelXD(powTrans(X))
#' # plots
#' curve(dnorm(x), col = "blue", xlim=c(-5, 5), ylim = c(0,0.7)) # true PDF
#' curve(dXD(x, xDens = xdens), add = TRUE, col = "red") # xDensity PDF
#' curve(dXD(x, xDens = xdens.trans), add = TRUE, col = "black") # xDensity PDF
#' legend("topleft", c("N(0, 1", "Chi-Sq(4)", "Power Trans"),
#'        pch = c(22,22,22,NA), pt.cex = 1.5,
#'        pt.bg = c("blue", "red", "black"))
#' @export
powTrans <- function(x, lambda = 0, alpha = 0, normalize = FALSE,
                      jacobian = FALSE) {
  if(lambda == 0) z <- log(x + alpha) else z <- ((x + alpha)^lambda - 1)/lambda
  if(normalize) {
    gm <- exp(mean(log(x)))
    if(lambda == 0) K <- gm else K <- 1/gm^(lambda-1)
  } else K <- 1
  ans <- z * K
  if(jacobian) ans <- list(z = ans, jacobian = (x + alpha)^(lambda - 1) * K)
  ans
}

#' Maximum likelihood estimation for the Generalized Box-Cox transformation.
#'
#' @param x Vector of random samples from target density.
#' @param alpha Optional value of the offset parameter.  \code{alpha = FALSE} sets \code{alpha = 1 - min(x)}, thereby guaranteeing that \code{z = x + alpha >= 1}.  This or any scalar value of \code{alpha} finds the conditional MLE as a function of \code{lambda} only.  \code{alpha = NA} finds the joint MLE over \code{(lambda,alpha)}.
#' @param interval Range of \code{lambda} values for one dimensional optimization.
#' @param ... Additional arguments to pass to \code{optimize} or \code{optim}, for 1- or 2-parameter optimization.
#' @details The likelihood for optimization is
#' \preformatted{L(lambda, alpha | x) = prod(dnorm(z(x | lambda, alpha)) *  |dz(x | lambda, alpha) / dx|)},
#' where \code{z(x | lambda, alpha)} is the Box-Cox transformation.
#' @return Vector of length two containing the fitted and/or known values of \code{(lambda, alpha)}.
powFit <- function(x, alpha = NA, interval = c(-5, 5), ...) {
 n <- length(x)
 mx <- min(x)
 fl <- function(lambda) {
   z <- powTrans(x = x, lambda = lambda, alpha = 0)
   s2 <- var(z)*(n-1)/n
   -n/2 * log(s2) + (lambda-1) * lx
 }
 fal <- function(theta) {
   if(theta[1] + mx <= 0) return(-Inf)
   z <- powTrans(x = x, lambda = theta[2], alpha = theta[1])
   s2 <- var(z)*(n-1)/n
   -n/2 * log(s2) + (theta[2]-1) * sum(log(x + theta[1]))
 }
 if(!is.na(alpha)) {
   # 1-parameter optimization
   if(is.logical(alpha) && !alpha) alpha <- 1 - mx
   x <- x + alpha
   lx <- sum(log(x))
   ans <- c(alpha = alpha,
            lambda = optimize(fl, interval = interval, maximum = TRUE)$maximum)
 } else {
   # 2-parameter optimization
   ans <- optim(par = c(1 - mx, 0), fn = fal, control = list(fnscale = -1, ...))
   if(ans$convergence != 0) stop("optim failed to converge.")
   ans <- ans$par
   names(ans) <- c("alpha", "lambda")
 }
 ans
}
