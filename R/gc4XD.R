#' Construct an \code{xDensity} representation of a combined Box-Cox/Gram-Charlier density approximation.
#'
#' @param x Vector of random samples from density to approximate.
#' @param lambda Exponent of Box-Cox transform.  If \code{NULL} it is estimated from \code{x}.  See Details.
#' @param alpha Offset of the Box-Cox transform.  Default is no offset.  See Details.
#' @param cmom Optional vector of first 4 central moments of \code{x}.  If \code{NULL} these are estimated from \code{x}.
#' @param trim Scalar between 0 and 1; removes the \code{trim} fraction of extreme values from \code{x} for estimation of \code{lambda} and \code{cmom}, which are very sensitive to outliers.  \code{trim = FALSE} does not trim any values.
#' @param n,from,to Specifies a grid of values on which to evaluate the density (see \code{\link[stats]{density}}).
#' @param ... Additional parameters to Box-Cox fitting function \code{\link{powFit}}.
#' @param mean,sd Optional mean and standard deviation for \code{xDensity} representation.
#' @details \code{x} is first standardized to \code{z = x/sd(x) - min(x/sd(x), from) + 1}, before the Box-Cox transform is applied.
#'
#' For details on the Box-Cox transformation and Gram-Charlier approximation, see \code{\link{powFit}} and \code{\link{dgc4}} respectively.
#' @return An \code{xDensity} object.
#' @seealso \code{\link{powFit}}, \code{\link{dgc4}}, \code{\link{xDensity}}.
#' @examples
#' # xDensity approximation to a noncentral-t distribution
#'
#' # true parameters
#' lambda <- rnorm(1) # noncentrality parameter
#' nu <- runif(1, 4, 6) # degrees of freedom
#'
#' # simulate data (note the small sample size)
#' x <- rt(500, df = nu, ncp = lambda)
#'
#' # xDensity approximation
#' xDensK <- kernelXD(x) # kernel smoothing
#' xDensG <- gc4XD(x) # gc4 approximation
#'
#' # true vs approximate PDFs
#' xlim <- qt(c(.005, .995), df = nu, ncp = lambda) # range for plot
#' curve(dt(x, df = nu, ncp = lambda),
#'       from = xlim[1], to = xlim[2], ylab = "Density")
#' curve(dXD(x, xDens = xDensK), add = TRUE, col = "red")
#' curve(dXD(x, xDens = xDensG), add = TRUE, col = "blue")
#' legend("topleft", legend = c("True PDF", "xDensity: kernel", "xDensity: gc4"),
#'        fill = c("black", "red", "blue"))
#' @export
gc4XD <- function(x, lambda = NULL, alpha = 0, cmom = NULL, trim = .01,
                  n = 512, from, to, mean, sd, ...) {
  weighted <- FALSE # depreciated
  x.sd <- stats::sd(x)
  x <- x/x.sd
  if(missing(from) | missing(to)) {
    bw <- bw.nrd0(x) * 3
    if(missing(from)) from <- min(x) - bw else from <- from/x.sd
    if(missing(to)) to <- max(x) + bw else to <- to/x.sd
  }
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
  matrixXD(XY = XY, mean = mean, sd = sd)
}
