#' Univariate 4th-order Gram-Charlier density approximation.
#'
#' @param x Vector of density quantiles.
#' @param cmom Vector of first 4 central moments of distribution.
#' @param rm.neg Logical; if \code{TRUE} clips the density at smallest calculated non-negative value.
#' @param log Logical; if \code{TRUE} returns approximating density on log scale.
#' @return Vector of density evaluations.
#' @examples
#' df <- 3
#' nsamples <- 100
#' X <- rchisq(nsamples, df=df) # iid samples from the Chi-Sq(df)
#'
#' # calculate central moments
#' cmom <- mean(X)
#' cmom <- c(cmom, mean((X-cmom)^2), mean((X-cmom)^3), mean((X-cmom)^4))
#'
#' curve(dgc4(x, cmom = cmom, rm.neg = TRUE), from = 0, to = 10)
#' @export
dgc4 <- function(x, cmom, rm.neg = TRUE, log = FALSE) {
  # get 3rd and 4th cumulant
  mu <- cmom[1]
  sigma <- sqrt(cmom[2])
  k3 <- cmom[3]
  k4 <- cmom[4] - 3*cmom[2]^2
  # standardize
  y <- (x - mu)/sigma
  # hermite polys
  H3 <- y^3 - 3*y
  H4 <- y^4 - 6*y^2 + 3
  sum.terms <- 1 + k3*H3/sigma^3/6 + k4*H4/sigma^4/24
  if(rm.neg) sum.terms <- pmax(sum.terms, min(sum.terms[sum.terms > 0]))
  dens <- dnorm(x, mean = mu, sd = sigma) * sum.terms
  if(log) dens <- log(dens)
  dens
}

# generalized box-cox transformation
# optionally include the jacobian vector, which converts transformed density
# values back to the original scale.
## pow.trans <- function(x, lambda = 0, alpha = 0, normalize = FALSE,
##                       jacobian = FALSE) {
##   if(lambda == 0) z <- log(x + alpha) else z <- ((x + alpha)^lambda - 1)/lambda
##   if(normalize) {
##     gm <- exp(mean(log(x)))
##     if(lambda == 0) K <- gm else K <- 1/gm^(lambda-1)
##   } else K <- 1
##   ans <- z * K
##   if(jacobian) ans <- list(z = ans, jacobian = (x + alpha)^(lambda - 1) * K)
##   ans
## }
