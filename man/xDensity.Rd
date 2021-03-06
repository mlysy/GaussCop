% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xDensity.R
\name{xDensity}
\alias{xDensity}
\alias{dXD}
\alias{pXD}
\alias{qXD}
\alias{rXD}
\title{Implementation of \code{d/p/q/r} functions for \code{xDensity} distributions.}
\usage{
dXD(x, xDens, log = FALSE)

pXD(q, xDens, lower.tail = TRUE, log.p = FALSE)

qXD(p, xDens, lower.tail = TRUE, log.p = FALSE)

rXD(n, xDens)
}
\arguments{
\item{x, q}{Vector of quantiles.}

\item{xDens}{Object of class \code{xDensity} representing the distribution.}

\item{log, log.p}{Logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.}

\item{lower.tail}{Logical; if \code{TRUE} (default), probabilities are \code{P[X <= x]} otherwise, \code{P[X > x]}.}

\item{p}{Vector of probabilities.}

\item{n}{Number of observations.}
}
\value{
For the underlying \code{xDensity} object, \code{dXD} gives the density, \code{pXD} gives the distribution function, \code{qXD} gives the quantile function and \code{rXD} generates \code{n} random values.
}
\description{
Implementation of \code{d/p/q/r} functions for \code{xDensity} distributions.
}
\details{
Extended density (or \code{xDensity}) objects provide a compact representation of arbitrary one-dimensional distributions defined on the real line.  That is, an \code{xDensity} object is a list with the following elements:
\itemize{
  \item \code{xrng}, \code{ndens}: range and number of gridpoints defining the main density region, i.e. \code{xseq = seq(xrng[1], xrng[2], len = xn)}.
  \item \code{ypdf}, \code{ylpdf}, \code{ycdf}: density, log-density, and cdf on the grid.
  \item \code{mean}, \code{sd}: mean and standard deviation of a Normal distribution to use outside the specified density range.
}
}
\examples{
# xDensity representation of a N(0,1) distribution

# construct the xDensity object using the known PDF dnorm
xseq <- seq(-4, 4, len = 500) # where to evaluate density
xDens <- matrixXD(cbind(xseq, dnorm(xseq)))

# check random sampling
x <- rXD(1e5, xDens = xDens)
hist(x, breaks = 100, freq = FALSE)
curve(dnorm, add = TRUE, col = "red")

# check PDF
x <- rnorm(5)
rbind(true = dnorm(x, log = TRUE), xDens = dXD(x, xDens, log = TRUE))

# check CDF
rbind(true = pnorm(x, log = TRUE), xDens = pXD(x, xDens, log = TRUE))

# check inverse-CDF
probs <- runif(5)
rbind(true = qnorm(probs), xDens = qXD(probs, xDens))
}
\seealso{
\code{\link{matrixXD}}, \code{\link{kernelXD}}, \code{\link{gc4XD}} for various \code{xDensity} object constructors.
}
