## ---- eval = FALSE, echo = FALSE-----------------------------------------
#  
#  # R code to render vignette
#  rmarkdown::render("xDensity-quicktut.Rmd")
#  

## ------------------------------------------------------------------------
require(GaussCop) # load package

# simulate data
nsamples <- 1e5
X <- log(rexp(nsamples)) # iid samples from the log-Exponential

xDens <- kernelXD(X) # basic xDensity constructor
names(xDens)

## ----xdens_test, echo = FALSE--------------------------------------------
# check that xDensity d/p/q/r functions match the true distribution,
# by plotting true PDF/CDF vs xDensity representation.
# xDens: xDensity object.
# dtrue/ptrue: single-variable function of the true pdf/cdf.
# dname: name of true distribution (for legend)
# dlpos, plpos: legend position for pdf/cdf plot.
xdens_test <- function(xDens, dtrue, ptrue, dname, dlpos, plpos) {
  op <- par(no.readonly = TRUE)
  par(mfrow = c(1,2), mar = c(4,4,2,.5)+.1)
  # pdf and sampling check
  nsamples <- 1e5
  Xsim <- rXD(nsamples, xDens = xDens) # xDensity: random sampling
  # extend xDens range to show endpoint matching
  xlim <- xDens$xrng
  xlim <- (xlim - mean(xlim)) * 1.10 + mean(xlim)
  hist(Xsim, breaks = 100, freq = FALSE, # xDensity random sample
       xlab = "x", main = "PDF", xlim = xlim)
  curve(dXD(x, xDens = xDens), add = TRUE, col = "red") # xDensity PDF
  curve(dtrue, add = TRUE, col = "blue") # true PDF
  abline(v = xDens$xrng, lty = 2) # grid endpoints
  legend(dlpos, parse(text = c(dname, "dXD", "rXD", "xrng")),
         pch = c(22,22,22,NA), pt.cex = 1.5,
         pt.bg = c("blue", "red", "white", "black"),
         lty = c(NA, NA, NA, 2))
  # cdf check
  curve(pXD(x, xDens), # xDensity CDF
        from = xlim[1], to = xlim[2], col = "red",
        xlab = "x", main = "CDF", ylab = "Cumulative Probability")
  curve(ptrue, add = TRUE, col = "blue") # true CDF
  abline(v = xDens$xrng, lty = 2) # grid endpoints
  legend(plpos, parse(text = c(dname, "pXD", "xrng")),
         pch = c(22,22,NA), pt.cex = 1.5,
         pt.bg = c("blue", "red", "black"),
         lty = c(NA, NA, 2))
  invisible(par(op))
}

## ---- eval = FALSE-------------------------------------------------------
#  dXD(x, xDens = xDens, log = TRUE) # log-PDF
#  pXD(x, xDens = xDens, lower.tail = FALSE) # survival function (1-CDF)
#  qXD(p = .5, xDens = xDens) # median
#  rXD(n, xDens = xDens) # random sample

## ---- fig.width = 10, fig.height = 4, out.width = "97%", fig.cap = "Graphical assessment of `xDensity` approximation to the log-Exponential distribution."----
# check that xDensity d/p/q/r functions match the true distribution,
# by plotting true PDF/CDF vs xDensity representation.

dlexp <- function(x) exp(-exp(x) + x) # true PDF
plexp <- function(x) pexp(exp(x)) # true CDF

xdens_test(xDens, # xDensity approximation 
           dtrue = dlexp, # true PDF 
           ptrue = plexp, # true CDF
           dname = expression(log(X)%~%"Expo"*(1)), # name of true distribution (for legend)
           dlpos = "topleft", plpos = "topleft") # legend position for pdf/cdf plot


## ---- fig.width = 10, fig.height = 4, out.width = "97%", fig.show = "hold", fig.cap = "Impact of sample size on quality of `kernelXD` approximation."----

# true distribution: chi^2_(4)
df_true <- 4
dtrue <- function(x) dchisq(x, df = df_true)
ptrue <- function(x) pchisq(x, df = df_true)
rtrue <- function(n) rchisq(n, df = df_true)
dname <- paste0("X%~%chi[(",df_true,")]^2")

# simulate data
nsamples <- 1e3
X <- rtrue(nsamples)
xDens <- kernelXD(X) # xDensity constructor

par(oma = c(0,2,0,0)) # put sample size in outer margin
xdens_test(xDens, dtrue = dtrue, ptrue = ptrue,
           dname = dname,
           dlpos = "topright", plpos = "bottomright")
mtext(text = paste0("M = ", nsamples),
      side = 2, line = .5, outer = TRUE)

# increase sample size
nsamples <- 1e4
X <- rtrue(nsamples)
xDens <- kernelXD(X)
xdens_test(xDens, dtrue = dtrue, ptrue = ptrue,
           dname = dname,
           dlpos = "topright", plpos = "bottomright")
mtext(text = paste0("M = ", nsamples),
      side = 2, line = .5, outer = TRUE)
par(oma = c(0,0,0,0)) # remove outer margin


## ---- echo = -1, fig.width = 10, fig.height = 4, out.width = "97%", fig.cap = "`gc4XD` approximation with $M = 1000$ samples."----
set.seed(1) # reproducible results
nsamples <- 1e3
X <- rtrue(nsamples)
xDens <- gc4XD(X) # xDensity constructor
xdens_test(xDens, dtrue = dtrue, ptrue = ptrue,
           dname = dname,
           dlpos = "topright", plpos = "bottomright")

## ---- fig.width = 10, fig.height = 4, out.width = "97%", fig.show = "hold", fig.cap = "`matrixXD` approximation for $t_{(\\nu)}$ distributions with $\\nu = 1,2$."----
# true distribution: t(2)
dtrue <- list(t2 = function(x) dt(x, df = 2),
              Cauchy = function(x) dt(x, df = 1))
ptrue <- list(t2 = function(x) pt(x, df = 2),
              Cauchy = function(x) pt(x, df = 1))
dname <- expression(X%~%t[(2)], X%~%"Cauchy")

# matrix density representation
xgrid <- seq(from=-15, to=15, len = 500) # grid
for(ii in 1:2) {
  ypdf <- dtrue[[ii]](xgrid) # true density
  xDens <- matrixXD(cbind(X = xgrid, Y = ypdf)) # xDensity constructor
  xdens_test(xDens, dtrue = dtrue[[ii]], ptrue = ptrue[[ii]],
             dname = dname[ii],
             dlpos = "topright", plpos = "bottomright")
}



## ---- fig.width = 10, fig.height = 2.5, out.width = "97%", fig.show = "hold", fig.cap = "Gaussian Copula approximation to the eigenvalues of $\\XX \\sim \\tx{Wishart}(\\II, 4)$, as fitted to $\\tx{chol}(\\XX)$."----
# plot true and approximate eigenvalue distributions
main <- parse(text = paste0("lambda[", 1:p, "]"))
par(mfrow = c(1,p), mar = c(2.5,2.5,2,.1)+.1, oma = c(0,2,0,0))
for(ii in 1:p) {
  dtrue <- density(Xeig[,ii])
  dapprox <- density(Xeig2[,ii])
  xlim <- range(dtrue$x, dapprox$x)
  ylim <- range(dtrue$y, dapprox$y)
  plot(dtrue$x, dtrue$y, type = "l", xlim = xlim, ylim = ylim,
       main = main[ii], xlab = "", ylab = "", col = "blue")
  lines(dapprox$x, dapprox$y, col = "red")
}
mtext(text = "Density",
      side = 2, line = .5, outer = TRUE)
legend("topright", legend = c("True PDF", "gCop Approx."),
       fill = c("blue", "red"))
par(oma = c(0,0,0,0)) # remove outer margin


## ---- ref.label = "xdens_test", eval = FALSE-----------------------------
#  # check that xDensity d/p/q/r functions match the true distribution,
#  # by plotting true PDF/CDF vs xDensity representation.
#  # xDens: xDensity object.
#  # dtrue/ptrue: single-variable function of the true pdf/cdf.
#  # dname: name of true distribution (for legend)
#  # dlpos, plpos: legend position for pdf/cdf plot.
#  xdens_test <- function(xDens, dtrue, ptrue, dname, dlpos, plpos) {
#    op <- par(no.readonly = TRUE)
#    par(mfrow = c(1,2), mar = c(4,4,2,.5)+.1)
#    # pdf and sampling check
#    nsamples <- 1e5
#    Xsim <- rXD(nsamples, xDens = xDens) # xDensity: random sampling
#    # extend xDens range to show endpoint matching
#    xlim <- xDens$xrng
#    xlim <- (xlim - mean(xlim)) * 1.10 + mean(xlim)
#    hist(Xsim, breaks = 100, freq = FALSE, # xDensity random sample
#         xlab = "x", main = "PDF", xlim = xlim)
#    curve(dXD(x, xDens = xDens), add = TRUE, col = "red") # xDensity PDF
#    curve(dtrue, add = TRUE, col = "blue") # true PDF
#    abline(v = xDens$xrng, lty = 2) # grid endpoints
#    legend(dlpos, parse(text = c(dname, "dXD", "rXD", "xrng")),
#           pch = c(22,22,22,NA), pt.cex = 1.5,
#           pt.bg = c("blue", "red", "white", "black"),
#           lty = c(NA, NA, NA, 2))
#    # cdf check
#    curve(pXD(x, xDens), # xDensity CDF
#          from = xlim[1], to = xlim[2], col = "red",
#          xlab = "x", main = "CDF", ylab = "Cumulative Probability")
#    curve(ptrue, add = TRUE, col = "blue") # true CDF
#    abline(v = xDens$xrng, lty = 2) # grid endpoints
#    legend(plpos, parse(text = c(dname, "pXD", "xrng")),
#           pch = c(22,22,NA), pt.cex = 1.5,
#           pt.bg = c("blue", "red", "black"),
#           lty = c(NA, NA, 2))
#    invisible(par(op))
#  }

## ---- eval = FALSE, include = FALSE--------------------------------------
#  # simulate data from unit Wishart distribution
#  nsamples <- 1e5
#  X <- replicate(nsamples, crossprod(matrix(rnorm(4),2,2)))
#  
#  # Cholesky decomposition
#  Xchol <- apply(X, 3, function(V) chol(V)[upper.tri(V,diag = TRUE)])
#  Xchol <- t(Xchol)
#  
#  # fit Gaussian Copula approximation
#  gCop <- gcopFit(X = Xchol, # iid sample from multivariate distribution
#                  fitXD = "kernel") # xDensity approximation to each margin
#  
#  # Compare approximation to true distribution of eigenvalues of X
#  # eigenvalues of 2x2 matrix
#  eig22 <- function(V) {
#    c1 <- V[1,1] + V[2,2] # trace
#    c2 <- V[1,1]*V[2,2] - V[1,2]^2 # determinant
#    c2 <- sqrt(.25 * c1^2 - c2)
#    c(.5 * c1 + c2, .5 * c1 - c2)
#  }
#  
#  Xeig <- t(apply(X, 3, eig22)) # true distribution
#  Xeig2 <- rgcop(nsamples, gCop = gCop) # Gaussian Copula approximation
#  Xeig2 <- t(apply(Xeig2, 1, function(U) {
#    # recover full matrix
#    V <- matrix(0,2,2)
#    V[upper.tri(V,diag=TRUE)] <- U
#    V <- crossprod(V)
#    eig22(V)
#  }))
#  main = expression(lambda[1],lambda[2])
#  
#  par(mfrow = c(1,2), mar = c(2.5,2.5,2,.1)+.1, oma = c(0,2,0,0))
#  for(ii in 1:2) {
#    dtrue <- density(Xeig[,ii])
#    dapprox <- density(Xeig2[,ii])
#    xlim <- range(dtrue$x, dapprox$x)
#    ylim <- range(dtrue$y, dapprox$y)
#    plot(dtrue$x, dtrue$y, type = "l", xlim = xlim, ylim = ylim,
#         main = main[ii], xlab = "", ylab = "")
#    lines(dapprox$x, dapprox$y, col = "red")
#  }
#  mtext(text = "Density",
#        side = 2, line = .5, outer = TRUE)
#  par(oma = c(0,0,0,0)) # remove outer margin

