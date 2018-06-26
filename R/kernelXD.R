#' Construct an \code{xDensity} representation of a kernel smoothing estimator.
#'
#' @param x Vector of samples from the underlying distribution.
#' @param n,from,to Optional arguments to \code{\link[stats]{density}} which are used to set the grid on which to evaluate the kernel estimator.
#' @param mean,sd Optional mean and standard deviation arguments for the extended density.  Default to the mean and standard deviation of \code{x}.
#' @param any0 Logical; if \code{FALSE} forces the support of the density to be the real line.
#' @param ... Additional arguments to \code{\link[stats]{density}}.
#' @return An \code{xDensity} object.
#' @examples
#' # xDensity approximation to a noncentral-t distribution
#'
#' # true parameters
#' lambda <- rnorm(1) # noncentrality parameter
#' nu <- runif(1, 4, 6) # degrees of freedom
#'
#' # simulate data
#' x <- rt(1e4, df = nu, ncp = lambda)
#'
#' # xDensity approximation
#' xDens <- kernelXD(x)
#'
#' # true vs approximate PDF
#' curve(dt(x, df = nu, ncp = lambda),
#'       from = min(x), to = max(x), ylab = "Density")
#' curve(dXD(x, xDens = xDens), add = TRUE, col = "red")
#' legend("topleft", legend = c("True PDF", "xDensity Approx."),
#'        fill = c("black", "red"))
#' @export
kernelXD <- function(x, n = 512, from, to, mean, sd, any0 = FALSE, ...) {
  dens <- density(x, n = n, from = from, to = to, ...)
  XY <- cbind(x = dens$x, y = dens$y)
  if(!any0) {
    # set zero's to smallest positive value
    XY[,2] <- pmax(XY[,2], min(XY[XY[,2] > 0,2]))
  }
  matrixXD(XY = XY, mean = mean, sd = sd)
}
