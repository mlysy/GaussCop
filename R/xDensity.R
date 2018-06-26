#' Implementation of \code{d/p/q/r} functions for \code{xDensity} distributions.
#'
#' @name xDensity
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param xDens Object of class \code{xDensity} representing the distribution.
#' @param log,log.p Logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \code{P[X <= x]} otherwise, \code{P[X > x]}.
#' @details Extended density (or \code{xDensity}) objects provide a compact representation of arbitrary one-dimensional distributions defined on the real line.  That is, an \code{xDensity} object is a list with the following elements:
#' \itemize{
#'   \item \code{xrng}, \code{ndens}: range and number of gridpoints defining the main density region, i.e. \code{xseq = seq(xrng[1], xrng[2], len = xn)}.
#'   \item \code{ypdf}, \code{ylpdf}, \code{ycdf}: density, log-density, and cdf on the grid.
#'   \item \code{mean}, \code{sd}: mean and standard deviation of a Normal distribution to use outside the specified density range.
#' }
#' @return For the underlying \code{xDensity} object, \code{dXD} gives the density, \code{pXD} gives the distribution function, \code{qXD} gives the quantile function and \code{rXD} generates \code{n} random values.
#' @seealso \code{\link{matrixXD}}, \code{\link{kernelXD}}, \code{link{gc4XD}} for various \code{xDensity} object constructors.
#' @examples
#' # xDensity representation of a N(0,1) distribution
#'
#' # construct the xDensity object using the known PDF dnorm
#' xseq <- seq(-4, 4, len = 500) # where to evaluate density
#' xDens <- matrixXD(cbind(xseq, dnorm(xseq)))
#'
#' # check random sampling
#' x <- rXD(1e5, xDens = xDens)
#' hist(x, breaks = 100, freq = FALSE)
#' curve(dnorm, add = TRUE, col = "red")
#'
#' # check PDF
#' x <- rnorm(5)
#' rbind(true = dnorm(x, log = TRUE), xDens = dXD(x, xDens, log = TRUE))
#'
#' # check CDF
#' rbind(true = pnorm(x, log = TRUE), xDens = pXD(x, xDens, log = TRUE))
#'
#' # check inverse-CDF
#' probs <- runif(5)
#' rbind(true = qnorm(probs), xDens = qXD(probs, xDens))
#' @rdname xDensity
#' @export
dXD <- function(x, xDens, log = FALSE) {
  if(class(xDens) != "xDensity") stop("xDens must be an xDensity object.")
  # extract components
  xrng <- xDens$xrng
  ndens <- xDens$ndens
  dx <- (xrng[2]-xrng[1])/ndens
  ylpdf <- xDens$ylpdf
  mean <- xDens$mean
  sd <- xDens$sd
  # calculate pdf
  ind <- (x - xrng[1])%/%dx + 1
  iind <- pmax(pmin(ind, ndens), 1)
  ans <- ifelse(1 <= ind & ind <= ndens,
                ylpdf[iind],
                dnorm(x, mean = mean, sd = sd, log = TRUE))
  if(!log) ans <- exp(ans)
  ans
}

#' @rdname xDensity
#' @export
pXD <- function(q, xDens, lower.tail = TRUE, log.p = FALSE) {
  if(class(xDens) != "xDensity") stop("xDens must be an xDensity object.")
  # extract components
  xrng <- xDens$xrng
  ndens <- xDens$ndens
  dx <- (xrng[2]-xrng[1])/ndens
  ypdf <- xDens$ypdf
  ycdf <- xDens$ycdf
  mean <- xDens$mean
  sd <- xDens$sd
  # calculate cdf
  ind <- (q - xrng[1])%/%dx + 1
  iind <- pmax(pmin(ind, ndens), 1)
  frac <- (q - xrng[1])%%dx
  if(lower.tail) {
    ans <- ifelse(1 <= ind & ind <= ndens,
                  ycdf[iind] + frac*ypdf[iind],
                  pnorm(q, mean = mean, sd = sd))
  } else {
    ans <- ifelse(1 <= ind & ind <= ndens,
                  1 - (ycdf[iind] + frac*ypdf[iind]),
                  pnorm(q, mean = mean, sd = sd, lower.tail = FALSE))
  }
  if(log.p) ans <- log(ans)
  ans
}

#' @rdname xDensity
#' @export
qXD <- function(p, xDens, lower.tail = TRUE, log.p = FALSE) {
  if(class(xDens) != "xDensity") stop("xDens must be an xDensity object.")
  # extract components
  ycdf <- xDens$ycdf
  xrng <- xDens$xrng
  ndens <- xDens$ndens
  dx <- (xrng[2]-xrng[1])/ndens
  xgrid <- seq(xrng[1]+.5*dx, xrng[2]-.5*dx, len = ndens)
  mean <- xDens$mean
  sd <- xDens$sd
  # calculate quantiles
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  tmp <- pmax(pmin(c(0, ycdf, 1), 1), 0)
  tmpi <- !duplicated(tmp)
  #ind <- as.numeric(cut(p, breaks = tmp[tmpi])) - 1
  ind <- findInterval(x = p, vec = tmp[tmpi]) - 1
  iind <- pmax(pmin(ind, sum(tmpi)-3), 1)
  Dy <- diff(ycdf[tmpi])[iind]
  PDy <- p - ycdf[tmpi][iind]
  ans <- ifelse(1 <= ind & ind <= sum(tmpi)-3,
                xgrid[tmpi][iind] + (PDy/Dy - 1/2)*dx,
                qnorm(p, mean = mean, sd = sd))
  #badp <- p < 0 | p > 1
  #if(any(badp)) {
  #  ans[badp] <- NaN
  #  warning("NaNs produced")
  #}
  ans
}

#' @rdname xDensity
#' @export
rXD <- function(n, xDens) {
  qXD(p = runif(n), xDens = xDens)
}

#--- examples ------------------------------------------------------------------

