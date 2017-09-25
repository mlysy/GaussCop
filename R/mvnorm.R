#--- internal multivariate normal functions ------------------------------------

.solveV <- function(V, x, ldV = FALSE) {
  C <- chol(V)
  if(missing(x)) x <- diag(nrow(V))
  ans <- backsolve(r = C, x = backsolve(r = C, x = x, transpose = TRUE))
  if(ldV) {
    ldV <- 2 * sum(log(diag(C)))
    ans <- list(y = ans, ldV = ldV)
  }
  ans
}

.lmvn <- function(x, sigma) {
  l2pi <- 1.837877066409345483560659472811 # log(2*pi)
  if(!is.matrix(x)) {
    tx <- as.matrix(x)
  } else {
    tx <- t(x)
  }
  iV <- .solveV(sigma, tx, ldV = TRUE)
  -.5 * (colSums(tx * iV$y) + iV$ldV + nrow(tx) * l2pi)
}

.rmvn <- function(n, sigma) {
  p <- ncol(sigma)
  matrix(rnorm(n*p), n, p) %*% chol(sigma)
}
