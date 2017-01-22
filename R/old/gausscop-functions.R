##########################################################################

# functions for gaussian copula estimation, smoothing, simulation, etc.

# copyright martin lysy, 2014
# mlysy@uwaterloo.ca

#############################################################

# gaussian copula functions
require(mvtnorm)

# parameter estimation of a Gaussian copula distribution.
cop.par <- function(X, dens.x, dens.y, Rho, mean, sd, from, to, n, ...,
                    zero.dens = FALSE, sd.infl = 1/5, debug = FALSE) {
  if(!missing(X)) {
    X <- as.matrix(X)
    d <- ncol(X)
    if(missing(mean)) mean <- colMeans(X)
    calc.mean <- FALSE
    if(missing(sd)) {
      sd <- 1
      sd <- apply(X, 2, sd)*sd.infl
    }
    calc.sd <- FALSE
    if(missing(Rho) & (missing(dens.x) | missing(dens.y))) {
      tmp <- apply(X, 2,
                   function(x) qqnorm(x, plot.it = FALSE)$x)
      Rho <- cor(tmp)
    }
    if(missing(dens.x) | missing(dens.y)) {
      if(!missing(from)) {
        from <- rep(from, len = d)
        to <- rep(to, len = d)
        if(missing(n)) n <- 512
        n <- rep(n, len = d)
        dens.x <- vector("list", d)
        names(dens.x) <- colnames(X)
        dens.y <- dens.x
        for(i in 1:d) {
          dens <- density(X[,i], n = n[i], from = from[i], to = to[i], ...)
          dens.x[[i]] <- dens$x
          dens.y[[i]] <- dens$y
        }
      } else {
        if(missing(n)) n <- 512
        n <- rep(n, len = d)
        dens.x <- vector("list", d)
        names(dens.x) <- colnames(X)
        dens.y <- dens.x
        for(i in 1:d) {
          dens <- density(X[,i], n = n[i], ...)
          dens.x[[i]] <- dens$x
          dens.y[[i]] <- dens$y
        }
      }
    }
  }
  if(debug) browser()
  # format arguments
  # dens.y
  if(is.matrix(dens.y)) {
    dens.y <- unlist(apply(dens.y, 2, function(x) list(x)), recursive = FALSE)
  } else if(is.numeric(dens.y)) dens.y <- list(dens.y)
  if(!is.list(dens.y)) stop("Incorrect specification for dens.y.")
  if(!zero.dens) dens.y <- sapply(dens.y, function(x) pmax(x, min(x[x > 0])),
                                  simplify = FALSE)
  ldens.y <- sapply(dens.y, log, simplify = FALSE)
  nrv <- length(dens.y)
  ndens <- sapply(dens.y, length)
  # dens.x
  if(is.matrix(dens.x)) {
    dens.x <- unlist(apply(dens.x, 2, function(x) list(x)), recursive = FALSE)
  } else if(is.numeric(dens.x)) dens.x <- list(dens.x)
  if(!is.list(dens.x)) stop("Incorrect specification for dens.x.")
  # space between dens evaluations
  dx <- sapply(dens.x, function(x) x[2]-x[1])
  # range of density: outside is normal approx
  rx <- sapply(dens.x, range)
  rx[1,] <- rx[1,] - dx/2
  rx[2,] <- rx[2,] + dx/2
  if(length(dens.x) != nrv)
    stop("nrv specified by dens.x and dens.y do not agree.")
  if(any(sapply(dens.x, length) != ndens))
    stop("ndens specified by dens.x and dens.y do not agree.")
  # mean and sd
  if(missing(mean)) {
    calc.mean <- TRUE
    mean <- rep(NA, nrv)
  } else {
    if(length(mean) != nrv) stop("mean has wrong length.")
    calc.mean <- FALSE
  }
  if(missing(sd)) {
    calc.sd <- TRUE
    sd <- rep(NA, nrv)
  } else {
    if(length(sd) != nrv) stop("sd has wrong length.")
    calc.sd <- FALSE
  }
  if(missing(Rho)) Z <- matrix(NA, nrow(X), ncol(X))
  # renormalize dens.y to account for tails; calculate CDF, mean, and sd
  Dens.y <- vector("list", nrv)
  names(Dens.y) <- names(dens.y)
  for(i in 1:nrv) {
    # calculate mean and sd
    if(calc.mean) mean[i] <- sum(dens.x[[i]]*dens.y[[i]])*dx[i]
    if(calc.sd)
      sd[i] <- sqrt(sum(dens.x[[i]]^2*dens.y[[i]])*dx[i] - mean[i]^2)*sd.infl
    # tail probability not contained in dens estimate
    tail.probs <- pnorm(abs(rx[,i] - mean[i])/sd[i], lower.tail = FALSE)
    dens.y[[i]] <- dens.y[[i]]/sum(dens.y[[i]]*dx[i])*(1-sum(tail.probs))
    Dens.y[[i]] <- cumsum(c(tail.probs[1], dens.y[[i]]*dx[i]))
    # calculate Rho
    if(missing(Rho)) {
      ind <- (X[,i] - rx[1,i]) %/% dx[i] + 1
      p <- (X[,i] - rx[1,i]) %% dx[i]
      iind <- pmax(pmin(ind, ndens[i]), 1)
      Z[,i] <- ifelse(1 <= ind & ind <= ndens[i],
                      qnorm(Dens.y[[i]][iind] + p*dens.y[[i]][iind]),
                      (X[,i] - mean[i])/sd[i])
    }
  }
  # Rho
  if(missing(Rho)) Rho <- cor(Z)
  Rho <- cov2cor(Rho)
  list(dens.x = dens.x, dx = dx, rx = rx,
       dens.y = dens.y, ldens.y = ldens.y, Dens.y = Dens.y,
       Rho = Rho, mean = mean, sd = sd)
}

# density of a Gaussian copula distribution
# decomp = TRUE returns the marginally Normalized data z, its multivariate
# log-density ldens.z, and the log-jacobian ljac.z, such that the total
# log-density is ldens.z + ljac.z
# max.Z: truncate normals to +/- max.Z standard deviations.
dcop <- function(x, par, dens.x, dens.y, Rho, mean, sd, log = FALSE,
                 decomp = FALSE, max.Z = 10, debug = FALSE) {
  if(missing(par)) par <- cop.par(dens.x = dens.x, dens.y = dens.y,
                                  Rho = Rho, mean = mean, sd = sd)
  dens.x <- par$dens.x
  dens.y <- par$dens.y
  ldens.y <- par$ldens.y
  Dens.y <- par$Dens.y
  Rho <- par$Rho
  mean <- par$mean
  sd <- par$sd
  dx <- par$dx
  rx <- par$rx
  nrv <- length(dens.x)
  ndens <- sapply(dens.x, length)
  # x
  if(!is.matrix(x)) {
    x <- as.matrix(x)
    if(nrv != 1) x <- t(x)
  }
  if(ncol(x) != nrv) stop("nrv specified by x and dens do not agree.")
  nx <- nrow(x)
  # pdf and cdf calculations
  f <- matrix(NA, nx, nrv)
  F <- f
  if(debug) browser()
  for(i in 1:nrv) {
    ind <- (x[,i] - rx[1,i])%/%dx[i] + 1
    iind <- pmax(pmin(ind, ndens[i]), 1)
    p <- (x[,i] - rx[1,i])%%dx[i]
    f[,i] <- ifelse(1 <= ind & ind <= ndens[i],
                    ldens.y[[i]][iind],
                    dnorm(x[,i], mean = mean[i], sd = sd[i], log = TRUE))
    F[,i] <- ifelse(1 <= ind & ind <= ndens[i],
                    Dens.y[[i]][iind] + p*dens.y[[i]][iind],
                    pnorm(x[,i], mean = mean[i], sd = sd[i]))
  }
  max.Z <- abs(max.Z)
  Z <- qnorm(F)
  Z <- pmax(pmin(Z, max.Z), -max.Z)
  ldens.z <- dmvnorm(Z, sigma = Rho, log = TRUE)
  jac.z <- rowSums(f) - dmvnorm(Z, sigma = diag(nrv), log = TRUE)
  if(decomp) return(list(Z = Z, ldens.z = ldens.z, jac.z = jac.z))
  ans <- ldens.z + jac.z
  if(!log) ans <- exp(ans)
  ans
}

# simulation of a Gaussian copula distribution.
# NOTE: problematic treatment of duplicate values in Dens.y.  Only seems to
# be a problem for highly skewed marginals, so far...
rcop <- function(n, par, dens.x, dens.y, Rho, mean, sd, debug = FALSE) {
  if(missing(par)) par <- cop.par(dens = dens.x, dens.y = dens.y,
                                  Rho = Rho, mean = mean, sd = sd)
  dens.x <- par$dens.x
  dens.y <- par$dens.y
  ldens.y <- par$ldens.y
  Dens.y <- par$Dens.y
  Rho <- par$Rho
  mean <- par$mean
  sd <- par$sd
  dx <- par$dx
  rx <- par$rx
  nrv <- length(dens.x)
  ndens <- sapply(dens.x, length)
  if(debug) browser()
  # Simulate Uniforms
  P <- pnorm(rmvnorm(n, sigma = Rho))
  X <- matrix(NA, n, nrv)
  nm <- which.min(sapply(par, function(l) is.null(names(l))))
  colnames(X) <- names(par[[nm]])
  for(i in 1:nrv) {
    tmp <- c(0, Dens.y[[i]], 1)
    tmpi <- !duplicated(tmp)
    ind <- as.numeric(cut(P[,i], breaks = tmp[tmpi])) - 1
    iind <- pmax(pmin(ind, sum(tmpi)-3), 1)
    Dy <- diff(Dens.y[[i]][tmpi])[iind]
    PDy <- P[,i] - Dens.y[[i]][tmpi][iind]
    X[,i] <- ifelse(1 <= ind & ind <= sum(tmpi)-3,
                    dens.x[[i]][tmpi][iind] + (PDy/Dy - 1/2)*dx[i],
                    qnorm(P[,i], mean = mean[i], sd = sd[i]))
  }
  X
}

#####################################################################3

# gram-charlier expansions

# gram charlier expansions.

# trimmed mean based not on quantile, but on contribution.
# i.e. eliminate the heaviest 1%.
# if weighted = TRUE, then elimate the top 1% of the weight contributing
# to the sum.
trimmed.mean <- function(x, trim = 0, weighted = FALSE, debug = FALSE) {
  if(trim == 0) return(mean(x))
  abs.x <- abs(x)
  if(debug) browser()
  if(!weighted) q <- quantile(abs.x, probs = 1-trim)
  else q <- wquantile(abs.x, weights = abs.x, probs = 1-trim)
  mean(x[abs.x <= q])
}

# generalized box-cox transformation
# optionally include the jacobian vector, which converts transformed density
# values back to the original scale.
pow.trans <- function(x, lambda = 0, alpha = 0, normalize = FALSE,
                      jacobian = FALSE, debug = FALSE) {
  if(lambda == 0) z <- log(x + alpha) else z <- ((x + alpha)^lambda - 1)/lambda
  if(debug) browser()
  if(normalize) {
    gm <- exp(mean(log(x)))
    if(lambda == 0) K <- gm else K <- 1/gm^(lambda-1)
  } else K <- 1
  ans <- z * K
  if(jacobian) ans <- list(z = ans, jacobian = (x + alpha)^(lambda - 1) * K)
  ans
}

# maximum likelihood estimate of lambda and alpha parameters
# alpha = FALSE fixes alpha = 1 - min(x), thereby guaranteeing that
# z = x + alpha >= 1 > 0.
# alpha = NA jointly maximizes lambda and alpha.
pow.mle <- function(x, alpha = NA, interval = c(-5, 5), ..., debug = FALSE) {
 n <- length(x)
 mx <- min(x)
 fl <- function(lambda) {
   z <- pow.trans(x = x, lambda = lambda, alpha = 0)
   s2 <- var(z)*(n-1)/n
   -n/2 * log(s2) + (lambda-1) * lx
 }
 fal <- function(theta) {
   if(theta[1] + mx <= 0) return(-Inf)
   z <- pow.trans(x = x, lambda = theta[2], alpha = theta[1])
   s2 <- var(z)*(n-1)/n
   -n/2 * log(s2) + (theta[2]-1) * sum(log(x + theta[1]))
 }
 if(debug) browser()
 if(!is.na(alpha)) {
   # 1-parameter optimization
   if(is.logical(alpha) && !alpha) alpha <- 1 - mx
   x <- x + alpha
   lx <- sum(log(x))
   ans <- c(alpha = alpha,
            lambda = optimize(fl, interval = interval, maximum = TRUE)$maximum)
 } else {
   # 2-parameter optimization
   ans <- optim(par = c(1 - mx, 0), fn = fal, control = list(fnscale = -1, ...))
   if(ans$convergence != 0) stop("optim failed to converge.")
   ans <- ans$par
   names(ans) <- c("alpha", "lambda")
 }
 ans
}

# univariate 4th order gram-charlier expansion
# rm.neg = TRUE puts the smallest non-negative value of the sum term in place
# of negative values.
dgc4 <- function(x, central.moments, log = FALSE, rm.neg = FALSE,
                 debug = FALSE) {
  # get 3rd and 4th cumulant
  mu <- central.moments[1]
  sigma <- sqrt(central.moments[2])
  k3 <- central.moments[3]
  k4 <- central.moments[4] - 3*central.moments[2]^2
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

# combined power transform + gram-charlier density estimated from a sample.
# x is first standardized to z = x/sd(x) - min(x/sd(x), from) + 1,
# before box-cox transform is applied.
# n, from, and to specify a grid of values to compute the density on.
# from and to default to the same values as the default kernel density method.
# lambda = NULL estimates the parameter by maximum likelihood, in which case
# ... are additional parameters to pass to pow.mle
# trim removes the outer fraction of values for computation of lambda and
# the moments, which are very sensitive to outliers.  trim = FALSE does not
# trim any values.
# output is a 2-column matrix of range and density values.
gc4.dens <- function(x, lambda = NULL, alpha = 0, central.moments = NULL,
                     trim = .01, weighted = FALSE,
                     n = 512, from, to, rm.neg = FALSE, ..., debug = FALSE) {
  x.sd <- sd(x)
  x <- x/x.sd
  if(missing(from) | missing(to)) {
    bw <- bw.nrd0(x) * 3
    if(missing(from)) from <- min(x) - bw else from <- from/x.sd
    if(missing(to)) to <- max(x) + bw else to <- to/x.sd
  }
  if(debug) browser()
  m <- min(x, from)
  z <- x - m + 1
  if(is.null(lambda)) lambda <- pow.mle(z, alpha = alpha, ...)["lambda"]
  z <- pow.trans(z, lambda = lambda, alpha = alpha)
  if(is.null(central.moments)) {
    cm <- mean(z)
    cm <- c(cm, trimmed.mean((z - cm)^2, trim = trim, weighted = weighted),
            trimmed.mean((z - cm)^3, trim = trim, weighted = weighted),
            trimmed.mean((z - cm)^4, trim = trim, weighted = weighted))
    central.moments <- cm
  }
  x <- seq(from, to, len = n)
  z <- pow.trans(x - m + 1, lambda = lambda, alpha = alpha, jacobian = TRUE)
  y <- dgc4(z$z, central.moments, rm.neg = rm.neg) * z$jacobian
  cbind(x = x*x.sd, y = y/x.sd)
}

# quantile function with weights
wquantile <- function(x, weights = NULL, probs = seq(0, 1, 0.25),
                      names = TRUE, debug = FALSE) {
  if(is.null(weights)) return(quantile(x, probs = probs, names = names))
  if(debug) browser()
  weights <- cumsum(weights[order(x)])
  weights <- weights/max(weights)
  x <- sort(x)
  # for some reason doesn't always work for p = 1
  q <- function(p) {
    if(p == 1) return(max(x))
    x[which.max(weights >= p)]
  }
  ans <- as.numeric(sapply(probs, q))
  if(names) names(ans) <- paste(probs*100, "%", sep = "")
  ans
}

###################################################

# depreciated
if(FALSE) {
pow.mle <- function(x, alpha = NA, interval = c(-5, 5),debug = FALSE) {
 n <- length(x)
 if(is.na(alpha)) alpha <- 1 - min(x)
 x <- x + alpha
 lx <- sum(log(x))
 f <- function(lambda) {
   z <- pow.trans(x, lambda = lambda, alpha = 0)
   s2 <- var(z)*(n-1)/n
   -n/2 * log(s2) + (lambda-1) * lx
 }
 if(debug) browser()
 c(alpha = alpha,
   lambda = optimize(f, interval = interval, maximum = TRUE)$maximum)
}
}
