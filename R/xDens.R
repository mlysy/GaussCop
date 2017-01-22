#' Extended-Density Objects.
#'
#' @name xDens
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param xDens object containing the extended-density (\code{xDens}) representation.
#' @param log,log.p logical; if TRUE, probabilities \code{p} are given as \code{log(p)}.
#' @param lower.tail logical; if TRUE (default), probabilities are \code{P[X <= x]} otherwise, \code{P[X > x]}.
#' @details An \code{xDens} object is a list with the following elements:
#' \itemize{
#'   \item \code{xrng}, \code{ndens}: range and number of gridpoints defining the main density region, i.e. \code{xseq = seq(xrng[1], xrng[2], len = xn)}.
#'   \item \code{ypdf}, \code{ylpdf}, \code{ycdf}: density, log-density, and cdf on the grid.
#'   \item \code{mean}, \code{sd}: mean and standard deviation of a normal to use outside the density range.
#' }

#' @rdname xDens
#' @export
dXD <- function(x, xDens, log = FALSE, debug = FALSE) {
  # extract components
  xrng <- xDens$xrng
  ndens <- xDens$ndens
  dx <- (xrng[2]-xrng[1])/ndens
  ylpdf <- xDens$ylpdf
  mean <- xDens$mean
  sd <- xDens$sd
  # calculate pdf
  if(debug) browser()
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
qXD <- function(p, xDens, lower.tail = TRUE, log.p = FALSE, debug = FALSE) {
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
  if(debug) browser()
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
