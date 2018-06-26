#' Construct an \code{xDensity} representation of a known distribution.
#'
#' @param XY 2-column matrix representation of a univariate density.  The first column is a grid of x values giving the centers of bins.   The second column are the corresponding densities.
#' @param mean,sd Optional mean and standard deviation of the Normal distribution to use outside the grid values in \code{XY}.  Defaults to grid-based approximation using \code{XY}.
#' @return An \code{xDensity} object.  See \code{\link{xDensity}}.
#' @examples
#' # xDensity approximation to a noncentral-t distribution
#'
#' # true parameters
#' lambda <- rnorm(1) # noncentrality parameter
#' nu <- runif(1, 4, 6) # degrees of freedom
#'
#' # xDensity encoding of known PDF
#' xlim <- qt(c(.005, .995), df = nu, ncp = lambda) # discretization range
#' xseq <- seq(xlim[1], xlim[2], len = 200) # where to evaluate density
#' xDens <- matrixXD(cbind(xseq, dt(xseq, df = nu, ncp = lambda)))
#'
#' # true vs approximate PDFs
#' curve(dt(x, df = nu, ncp = lambda),
#'       from = xlim[1], to = xlim[2], ylab = "Density")
#' curve(dXD(x, xDens = xDens), add = TRUE, col = "red")
#' legend("topleft", legend = c("True PDF", "xDensity Approx."),
#'        fill = c("black", "red"))
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
  class(xDens) <- "xDensity"
  xDens
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
