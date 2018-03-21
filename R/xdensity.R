#' @title Extend a matrix representation of a density to an \code{xDens} object.
#' @description Constructs an \code{xDens} object using a pre-specified matrix of grid and density values.
#' @param XY 2-column matrix representation of a univariate density.  The first column is a grid of x values giving the centers of bins.   The second column are the corresponding densities.
#' @param mean,sd Optional mean and standard deviation of the Normal distribution to use outside the grid values in \code{XY}.  Defaults to grid-based approximation using \code{XY}.
#' @return An \code{xDens} object.  See \code{\link{xDens}}.
#' @examples
#' # simulate data from Cauchy(0, 1) distribution
#' nsamples = 1e4
#' X = seq(from=-5, to=4, by=0.1) # grid
#' Y = dcauchy(X)
#' xDens = matrixXD(cbind(X, Y))
#'par(mfrow = c(1,2), mar = c(4,4,2,.5)+.1)
#' 
#' # pdf and sampling check
#' Xsim <- rXD(nsamples, xDens = xDens) # xDensity: random sampling
#' hist(Xsim, breaks = 100, freq = FALSE,
#'      xlab = "x", main = "PDF")
#' curve(dXD(x, xDens = xDens), add = TRUE, col = "red") # xDensity PDF
#' curve(dcauchy(x), add = TRUE, col = "blue") # true PDF
#' abline(v = xDens$xrng, lty = 2) # grid endpoints
#' legend("topleft", c("rXD", "dXD", "Cauchy(0,1)", "xRange"),
#'        pch = c(22,22,22,NA), pt.cex = 1,
#'        pt.bg = c("white", "red", "blue", "black"),
#'        lty = c(NA, NA, NA, 2))
#' 
#' # cdf check
#' curve(pXD(x, xDens), # xDens CDF
#'       from = min(Xsim), to = max(Xsim), col = "red",
#'       xlab = "x", main = "CDF", ylab = "Cumulative Probability")
#' curve(pcauchy(x), add = TRUE, col = "blue") # true CDF
#' abline(v = xDens$xrng, lty = 2) # grid endpoints
#' legend("bottomright", c("pXD", "Cauchy(0,1)", "xRange"),
#'        pch = c(22,22,NA), pt.cex = 1.5,
#'        pt.bg = c("red", "blue", "black"),
#'        lty = c(NA, NA, 2))
#' @export
matrixXD <- function(XY, mean, sd) {
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
  if(missing(mean) || missing(sd)) {
    summary.stats <- .getNorm(XY)
    if(missing(mean)) mean <- summary.stats[1]
    if(missing(sd)) sd <- summary.stats[2]
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

#' @title Extended Kernel Density Approximation Constructor.
#' @description Constructs an \code{xDens} object using kernel density approximation.
#' @param x Vector of samples from the underlying distribution.
#' @param n,from,to Optional arguments to \code{\link{stat::density}} which are used to set the grid on which to evaluate the kernel estimator.
#' @param mean,sd Optional mean and standard deviation arguments for the extended density.  Default to the mean and standard deviation of \code{x}.
#' @param any0 Logical; if FALSE forces the support of the density to be the real line.
#' @param ... additional arguments to \code{\link{stat::density}}.
#' @return an \code{xDens} object.
#' @examples 
#' df = 4 # chisq degrees of freedom
#' nsamples = 5e3
#' X <- rchisq(nsamples, df=df) # iid samples from the Chi-Sq(3)
#' xDens <- kernelXD(X) # Kernel density estimation constructor
#' par(mfrow = c(1,2), mar = c(4,4,2,.5)+.1)
#' 
#' # pdf and sampling check
#' Xsim <- rXD(nsamples, xDens = xDens) # xDensity: random sampling
#' hist(Xsim, breaks = 100, freq = FALSE,
#'      xlab = "x", main = "PDF", ylim=c(0, 0.2))
#' curve(dXD(x, xDens = xDens), add = TRUE, col = "red") # xDensity PDF
#' curve(dchisq(x, df=df), add = TRUE, col = "blue") # true PDF
#' abline(v = xDens$xrng, lty = 2) # grid endpoints
#' legend("topright", c("rXD", "dXD", "dchisq(4)", "xRange"),
#'        pch = c(22,22,22,NA), pt.cex = 1.5,
#'        pt.bg = c("white", "red", "blue", "black"),
#'        lty = c(NA, NA, NA, 2))
#' 
#' # cdf check
#' curve(pXD(x, xDens), # xDensity CDF
#'       from = min(Xsim), to = max(Xsim), col = "red",
#'       xlab = "x", main = "CDF", ylab = "Cumulative Probability")
#' curve(pchisq(x, df=df), add = TRUE, col = "blue") # true CDF
#' abline(v = xDens$xrng, lty = 2) # grid endpoints
#' legend("bottomright", c("pXD", "chisq(4)", "xRange"),
#'        pch = c(22,22,NA), pt.cex = 1.5,
#'        pt.bg = c("red", "blue", "black"),
#'        lty = c(NA, NA, 2))
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

#--- Gaussian Endpoint Matching ------------------------------------------------

## @param XY 2-column matrix representation of a univariate density.  The first column is a grid of x values giving the centers of bins. The second column are the corresponding densities.
## @return A vector of length 2 with the first element being the variance and second being the mean
## DON'T @export
.getNorm <- function(XY) {
  xgrid <- XY[,1]
  ypdf <- XY[,2]
  dx <- xgrid[2] - xgrid[1]
  # match on first/last non-zero elements of ypdf
  ig0 <- ypdf > 0
  ng0 <- sum(ig0)
  ylpdf <- log(ypdf[ig0][c(1, ng0)])
  xrng <- xgrid[ig0][c(1, ng0)]
  # endpoint matching constants
  a <- .5 * (xrng[1] + xrng[2])
  B <- (ylpdf[1] - ylpdf[2])/(xrng[1] - xrng[2])
  A <- (xrng[1] - a)
  C <- 2*B*(A) - 2*ylpdf[1] - log(2*pi)
  D <- 1 + log((A^2)/2)
  xa2 <- xrng[1] - a
  w <- 2*B*xa2 - 2*ylpdf[1] - log(2*pi)
  B2 <- B^2
  xa2 <- xa2^2
  ## return(c(a = a, B = B, A = A, C = C, D = D))
  # root-finding function
  g <- function (tau) {
    ## x <- xrng[1]
    ## y <- ylpdf[1]
    ## xa <- x - a
    ## w <- 2*B*xa - 2*y - log(2*pi)
    return(log(tau) + xa2/tau + B2*tau - w)
  }
  # root bounds
  U <- (A^2)/(2*(C-D)) # upper
  # lower bound
  if(B == 0) {
    L <- A^2
  } else {
    L <- 1 + 4*(A*B)^2
    L <- if(L < 0) 0 else (-1 + sqrt(L)) / (2*B^2)
  }
  ## return(c(L = L, U = U))
  if(L == 0 || g(L) > 0) {
    # check if solution exists: is this necessary?
    warning("Endpoint matching failed.  mean/sd set to empirical estimates.")
    mu <- sum(xgrid*ypdf)*dx # mean calculated from XY
    sig2 <- sum((xgrid - mu)^2*ypdf)*dx # variance calculated from XY
  } else {
    sig2 <- uniroot(g, interval = c(L, U))$root
    mu <- a + B*sig2[1]
  }
  return(c(mu, sqrt(sig2)))
}

#--- test code -----------------------------------------------------------------

## XY <- density(rchisq(1e5, df = 4))
## XY <- cbind(X = XY$x, Y = XY$y)

## getNorm(XY)
## getNorm_test(XY)

## getNorm_test <- function (XY) {
##   xgrid <- XY[,1]
##   dx <- diff(xgrid)[1]
##   ndens <- length(xgrid)
##   xrng <- range(xgrid) + c(-.5, .5) * dx
##   ypdf = XY[,2]
##   # check to see if endpoints are 0 before logging
##   if (ypdf[1] > 0 && ypdf[ndens] > 0) {
##     ylpdf <- log(c(ypdf[1], ypdf[ndens]))
##   }
##   else { # one of the endpoint densities = 0
##     warning("Grid endpoint density is zero. Density mean and standard deviation used in normally distributed tails")
##     mu.temp = sum(xgrid*ypdf)*dx # mean calculated from XY
##     tau.temp = sum((xgrid - mu.temp)^2*ypdf)*dx # variance calculated from XY
##     return(c(tau.temp, mu.temp))
##   }

##   a = mean(c(xrng[1], xrng[2]))
##   B = (ylpdf[1] - ylpdf[2])/(xrng[1] - xrng[2])

##   # function to find root of:
##   g <- function (tau) {
##     x = xrng[1]
##     y = ylpdf[1]
##     xa = x - a
##     w = 2*B*xa - 2*y - log(2*pi)
##     return (log(tau) + (xa^2)/tau + ((B^2)*tau) - w)
##   }

##   A = (xrng[1] - a)
##   C = 2*B*(A) - 2*ylpdf[1] - log(2*pi)
##   D = 1 + log((A^2)/2)
##   #return(c(a = a, B = B, A = A, C = C, D = D))


##   upper.bound = (A^2)/(2*(C-D))
##   lower.bound = (-1 + sqrt(1 + 4*(A*B)^2)) / (2*B^2)
##   ## return(c(L = lower.bound, U = upper.bound))

##   # plot g() function for testing
##   # curve(g(x), from=upper.bound-5, to = lower.bound + 50)
##   # abline(v=lower.bound, lty=2)
##   # abline(v=upper.bound, lty=2)
##   # abline(h=g(upper.bound), lty=2) # lower bound captured
##   # abline(h=g(lower.bound), lty=2) # lower bound captured

##   no.solution = FALSE
##   if (lower.bound == 0) { no.solution = TRUE } # ylpdf endpoint densities are too close causing B = 0
##   else {if(g(lower.bound) > 0) {no.solution = TRUE}} # no change in sign if g(lower.bound) > 0

##   if (no.solution) {  # no solution
##     warning("Density mean and standard deviation used in normally distributed tails")
##     mu.temp = sum(xgrid*ypdf)*dx # mean calculated from XY
##     tau.temp = sum((xgrid - mu.temp)^2*ypdf)*dx # variance calculated from XY
##     return(c(mu.temp, sqrt(tau.temp)))
##   }

##   sigma2.mu = uniroot(g, interval=c(lower.bound, upper.bound))$root
##   sigma2.mu <- c(a + B*sigma2.mu[1], sqrt(sigma2.mu))
##   ## sigma2.mu[2] = a + B*sigma2.mu[1]
##   return (sigma2.mu)
## }
