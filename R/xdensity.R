#' Extend a Matrix Representation of a Density to an \code{xDens} Object
#'
#' @param XY 2-column matrix representation of a univariate density.  The first column is a grid of x values giving the centers of bins.   The second column are the corresponding densities.
#' @param mean,sd mean and standard deviation of the normal to use outside the grid values in \code{XY}.
#' @return An \code{xDens} object.  See \code{\link{xDens}}.
#' @export
xdensity <- function(XY, mean, sd) {
  #sd.infl <- 1/5 # should eventually remove this...
  sd.infl <- 1
  xgrid <- XY[,1]
  # check for regular grid
  dx <- diff(xgrid)
  if(stats::sd(dx)/base::mean(dx) > 1e-6) {
    stop("First column of XY must be a regular grid.")
  }
  ndens <- length(xgrid)
  dx <- dx[1]
  xrng <- range(xgrid) + c(-.5, .5) * dx
  ypdf <- XY[,2]
  # mean and sd
  if(missing(mean)) {
    mean <- sum(xgrid*ypdf)*dx
  }
  if(missing(sd)) {
    sd <- sqrt(sum((xgrid - mean)^2*ypdf)*dx) * sd.infl
  }
  # tail probability not contained in dens estimate
  tprob <- pnorm(abs(xrng - mean)/sd, lower.tail = FALSE)
  # renormalize ypdf
  ypdf <- ypdf/sum(ypdf*dx)*(1-sum(tprob))
  ycdf <- cumsum(c(tprob[1], ypdf*dx))
  ylpdf <- log(ypdf)
  xDens <- list(ndens = ndens, xrng = xrng, ypdf = ypdf, ylpdf = ylpdf,
                ycdf = ycdf, mean = mean, sd = sd)
  class(xDens) <- "xDens"
  xDens
}

#' Extended Kernel Density Approximation
#'
#' @param x vector of samples from the underlying distribution.
#' @param n,from,to optional arguments to \code{\link{stat::density}} which are used to set the grid on which to evaluate the kernel estimator.
#' @param mean,sd optional mean and standard deviation arguments for the extended density.  Default to the mean and standard deviation of \code{x}.
#' @param any0 logical; if FALSE forces the support of the density to be the real line.
#' @param ... additional arguments to \code{\link{stat::density}}.
#' @return an \code{xDens} object.
#' @export
kernelXD <- function(x, n = 512, from, to, mean, sd, any0 = FALSE, ...) {
  dens <- density(x, n = n, from = from, to = to, ...)
  XY <- cbind(x = dens$x, y = dens$y)
  if(!any0) {
    # set zero's to smallest positive value
    XY[,2] <- pmax(XY[,2], min(XY[XY[,2] > 0,2]))
  }
  xdensity(XY = XY, mean = mean, sd = sd)
}

#' Normal Extended Density.
#'
#' @param x vector of samples to which to fit the normal distribution.
#' @param n,from,to optional arguments to set up the grid on which to store the density estimate.
#' @param mean,sd optional mean and standard deviation arguments for the extended density.  Default to the mean and standard deviation of \code{x}.
#' @return an \code{xDens} object.
#' @export
normalXD <- function(x, n = 512, from, to, mean, sd) {
  if(missing(from) | missing(to)) {
    bw <- bw.nrd0(x) * 3
    if(missing(from)) from <- min(x) - bw
    if(missing(to)) to <- max(x) + bw
  }
  if(missing(mean)) mean <- base::mean(x)
  if(missing(sd)) sd <- stats::sd(x)
  xgrid <- seq(from, to, len = n)
  xdensity(XY = cbind(x = xgrid,
             y = dnorm(x = xgrid, mean = mean, sd = sd)),
           mean = mean, sd = sd)
}
