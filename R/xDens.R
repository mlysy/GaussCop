#' @title Calculation of \code{d/p/q/r} functions for density approximations.
#' @description Extended Density, \code{xDens}, objects provide a unified framework for the calculation of R's standard \code{d/q/p/r} functions
#' for unconstrained density approximations of a finite set of values. See the `constructors` section of the vignette 
#' for more information on the construction of \code{xDens} objects.
#' @name xDens
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param xDens Object containing the extended-density (\code{xDens}) representation.
#' @param log,log.p Logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \code{P[X <= x]} otherwise, \code{P[X > x]}.
#' @details An \code{xDens} object is a list with the following elements:
#' \itemize{
#'   \item \code{xrng}, \code{ndens}: range and number of gridpoints defining the main density region, i.e. \code{xseq = seq(xrng[1], xrng[2], len = xn)}.
#'   \item \code{ypdf}, \code{ylpdf}, \code{ycdf}: density, log-density, and cdf on the grid.
#'   \item \code{mean}, \code{sd}: mean and standard deviation of a Normal distribution to use outside the specified density range.
#' }
#' @rdname xDens
#' @return For the underlying \code{xDens} object, \code{dXD} gives the density, \code{pXD} gives the distribution function, 
#' \code{qXD} gives the quantile function and \code{rXD} generates \code{n} random values. 
#' @examples 
#' nsamples <- 1e5
#' X <- rnorm(nsamples) # generate data
#' xDens <- kernelXD(X) # Use kernel constructor to construct xDens object
#' 
#' # pdf and sampling check
#' Xsim <- rXD(nsamples, xDens = xDens) # xDensity sampling using rXD
#' hist(Xsim, breaks = 100, freq = FALSE,
#'      xlab = "x", main = "PDF")
#' curve(dXD(x, xDens = xDens), add = TRUE, col = "red") # xDensity PDF
#' curve(dnorm(x), add = TRUE, col = "blue") # true PDF
#' abline(v = xDens$xrng, lty = 2) # grid endpoints
#' legend("topright", c("rXD", "dXD", "N(0,1)", "xRange"),
#'        pch = c(22,22,22,NA), pt.cex = 1.5,
#'        pt.bg = c("white", "red", "blue", "black"),
#'        lty = c(NA, NA, NA, 2))
#' 
#' # cdf check
#' curve(pXD(x, xDens), # xDensity CDF
#'       from = min(Xsim), to = max(Xsim), col = "red",
#'       xlab = "x", main = "CDF", ylab = "Cumulative Probability")
#' curve(pnorm(x), add = TRUE, col = "blue") # true CDF
#' abline(v = xDens$xrng, lty = 2) # grid endpoints
#' legend("bottomright", c("pXD", "N(0,1)", "xRange"),
#'        pch = c(22,22,NA), pt.cex = 1.5,
#'        pt.bg = c("red", "blue", "black"),
#'        lty = c(NA, NA, 2))
#'
#' # quantile check
#' curve(qXD(x, xDens), # xDensity quantile function
#'       from = 0, to = 1, col = "red",
#'       xlab = "Quantiles", main = "Quantile Plot", ylab = "x")
#' curve(qnorm(x), add = TRUE, col = "blue") # true CDF
#' legend("bottomright", c("qXD", "N(0, 1)"),
#'        pch = c(22,22,NA), pt.cex = 1.5,
#'        pt.bg = c("red", "blue"),
#'        lty = c(NA, NA, 2))
#' @export
dXD <- function(x, xDens, log = FALSE) {
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

#' @rdname xDens
#' @export
pXD <- function(q, xDens, lower.tail = TRUE, log.p = FALSE) {
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

#' @rdname xDens
#' @export
qXD <- function(p, xDens, lower.tail = TRUE, log.p = FALSE) {
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

#' @rdname xDens
#' @export
rXD <- function(n, xDens) {
  qXD(p = runif(n), xDens = xDens)
}

