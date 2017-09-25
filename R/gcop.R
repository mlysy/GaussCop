#' The Gaussian Copula Distribution.
#'
#' @name gcop
#' @param X \code{n x p} values at which to evaluate the \code{p}-dimensional density.
#' @param n number of random samples to draw.
#' @param gCop the \code{gcop} model specification.
#' @param decomp logical; if TRUE returns the normalized residuals \code{Z}, their log-density \code{zlpdf}, and the log-jacobian \code{zljac}, such that the total log-density is \code{xldens = zldens + zljac}.  See details.
#' @details The density of Gaussian Copula distribution is
#' \deqn{
#' g(x) = \frac{\psi(z \mid \Rho) \prod_{i=1}^d f_i(x_i)}{\prod_{i=1}^d\phi(z_i)},
#' }{
#' g(x) = \psi(z | R) \prod_{i=1}^d f_i(x_i)/\phi(z_i),
#' }
#' \deqn{
#' z_i = \Phi^{-1}(F_i(x_i)),
#' }{
#' z_i = \Phi^(-1)(F_i(x_i)),
#' }
#' where \eqn{\psi(z \mid \Rho)}{\psi(z | R)} is the PDF of a multivariate normal with mean 0 and variance \eqn{\Rho}{R}, \eqn{f_i(x_i)} and \eqn{F_i(x_i)} are the marginal PDF and CDF of variable \eqn{i}, and \eqn{\phi(z)} and \eqn{\Phi(z)} are the PDF and CDF of a standard normal.
#' @return A vector of densities or matrix of random samples.

#' @rdname gcop
#' @export
dgcop <- function(X, gCop, log = FALSE, decomp = FALSE) {
  maxZ <- 10 # truncate normals to +/- maxZ standard deviations.
  if(class(gCop) != "gcop") stop("gCop must be a gcop object.")
  # format X
  nrv <- length(gCop$XDens)
  if(!is.matrix(X)) {
    X <- as.matrix(X)
    if(nrv != 1) X <- t(X)
  }
  if(ncol(X) != nrv) {
    stop("number of variables in X and gCop don't agree.")
  }
  # lpdf and cdf evaluations
  lf <- matrix(NA, nrow(X), nrv)
  FF <- lf
  for(ii in 1:nrv) {
    lf[,ii] <- dXD(X[,ii], gCop$XDens[[ii]], log = TRUE)
    FF[,ii] <- pXD(X[,ii], gCop$XDens[[ii]])
  }
  # normalized quantiles
  Z <- qnorm(FF)
  Z <- pmax(pmin(Z, maxZ), -maxZ) # truncate to reasonable range
  zlpdf <- .lmvn(Z, sigma = gCop$Rho)
  zljac <- rowSums(lf - dnorm(Z, log = TRUE))
  # zljac <- rowSums(lf) - dmvnorm(Z, sigma = diag(nrv), log = TRUE)
  if(decomp) {
    return(list(Z = Z, zlpdf = zlpdf, zljac = zljac))
  }
  ans <- zlpdf + zljac
  if(!log) ans <- exp(ans)
  ans
}

#' @rdname gcop
#' @export
rgcop <- function(n, gCop) {
  if(class(gCop) != "gcop") stop("gCop must be a gcop object.")
  # simulate correlated uniforms
  nrv <- length(gCop$XDens)
  U <- pnorm(.rmvn(n, sigma = gCop$Rho))
  # invert cdfs to obtain marginals
  X <- matrix(NA, n, nrv)
  colnames(X) <- names(gCop$XDens)
  for(ii in 1:nrv) {
    X[,ii] <- qXD(p = U[,ii], xDens = gCop$XDens[[ii]])
  }
  X
}
