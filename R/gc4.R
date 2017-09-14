#' Univariate 4th Order Gram-Charlier Density
#'
#' @param x vector of density quantiles.
#' @param cmom vector of first 4 central moments of distribution
#' @param rm.neg logical; if TRUE clips the density at smallest calculated non-negative value.
#' @param log logical; if TRUE returns approximating density on log scale.
#' @return Vector of pdf or log(pdf).
#' @export
dgc4 <- function(x, cmom, rm.neg = TRUE, log = FALSE,
                 debug = FALSE) {
  # get 3rd and 4th cumulant
  mu <- cmom[1]
  sigma <- sqrt(cmom[2])
  k3 <- cmom[3]
  k4 <- cmom[4] - 3*cmom[2]^2
  # standardize
  y <- (x - mu)/sigma
  # hermite polys
  H3 <- y^3 - 3*y
  H4 <- y^4 - 6*y^2 + 3
  sum.terms <- 1 + k3*H3/sigma^3/6 + k4*H4/sigma^4/24
  if(rm.neg) sum.terms <- pmax(sum.terms, min(sum.terms[sum.terms > 0]))
  dens <- dnorm(x, mean = mu, sd = sigma) * sum.terms
  if(log) dens <- log(dens)
  dens
}

#' Combined Box-Cox Transformation + Gram-Charlier Density Approximation.
#'
#' @param x Sample from density to approximate.
#' @param lambda Exponent of Box-Cox transform.  If NULL it is estimated from \code{x}.  See details.
#' @param alpha Offset of the Box-Cox transform.  Default is no offset.  See details.
#' @param cmom Optional vector of first 4 central moments.  If NULL these are estimated from \code{x}.
#' @param trim Scalar between 0 and 1; removes the \code{trim} fraction of extreme values from \code{x} for estimation of \code{lambda} and \code{cmom}, which are very sensitive to outliers.  \code{trim = FALSE} does not trim any values.
#' @param n,from,to specifies a grid of values on which to evaluate the density.
#' @param ... additional parameters to Box-Cox fitting function \code{\link{powFit}}.
#' @param mean,sd mean and standard deviation for extended density.
#' @details \code{x} is first standardized to \code{z = x/sd(x) - min(x/sd(x), from) + 1}, before the Box-Cox transform is applied.
#'
#' \code{lambda} can be estimated from the data, but for stability \code{alpha} must be provided.
#' @return An \code{xDens} object.
#' @export
gc4XD <- function(x, lambda = NULL, alpha = 0, cmom = NULL, trim = .01,
                  n = 512, from, to, mean, sd, ..., debug = FALSE) {
  weighted <- FALSE # depreciated
  x.sd <- stats::sd(x)
  x <- x/x.sd
  if(missing(from) | missing(to)) {
    bw <- bw.nrd0(x) * 3
    if(missing(from)) from <- min(x) - bw else from <- from/x.sd
    if(missing(to)) to <- max(x) + bw else to <- to/x.sd
  }
  if(debug) browser()
  m <- min(x, from)
  z <- x - m + 1
  if(is.null(lambda)) lambda <- powFit(z, alpha = alpha, ...)["lambda"]
  z <- powTrans(z, lambda = lambda, alpha = alpha)
  if(is.null(cmom)) {
    cm <- base::mean(z)
    cm <- c(cm, trimmed.mean((z - cm)^2, trim = trim,
                             weighted = weighted),
            trimmed.mean((z - cm)^3, trim = trim, weighted = weighted),
            trimmed.mean((z - cm)^4, trim = trim, weighted = weighted))
    cmom <- cm
  }
  x <- seq(from, to, len = n)
  z <- powTrans(x - m + 1, lambda = lambda, alpha = alpha,
                jacobian = TRUE)
  y <- dgc4(z$z, cmom = cmom, rm.neg = TRUE) * z$jacobian
  XY <- cbind(x = x*x.sd, y = y/x.sd)
  xdensity(XY = XY, mean = mean, sd = sd)
}

# generalized box-cox transformation
# optionally include the jacobian vector, which converts transformed density
# values back to the original scale.
## pow.trans <- function(x, lambda = 0, alpha = 0, normalize = FALSE,
##                       jacobian = FALSE, debug = FALSE) {
##   if(lambda == 0) z <- log(x + alpha) else z <- ((x + alpha)^lambda - 1)/lambda
##   if(debug) browser()
##   if(normalize) {
##     gm <- exp(mean(log(x)))
##     if(lambda == 0) K <- gm else K <- 1/gm^(lambda-1)
##   } else K <- 1
##   ans <- z * K
##   if(jacobian) ans <- list(z = ans, jacobian = (x + alpha)^(lambda - 1) * K)
##   ans
## }
