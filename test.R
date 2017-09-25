#--- package testing ------------------------------------------------------------

require(GaussCop)
require(mvtnorm)

# simulate some data from a Gaussian copula

# marginal models
#rf1 <- function(n) log(rgamma(n, shape = 1))
#rf2 <- function(n) rt(n, df = 10)
#rf3 <- function(n) rbeta(n, shape1 = 1, shape2 = 4)
qf1 <- function(p) log(qgamma(p = p, shape = 1))
qf2 <- function(p) qt(p, df = 10)
qf3 <- function(p) qbeta(p, shape1 = 1, shape2 = 4)

# correlation
Rho <- cov2cor(crossprod(matrix(rnorm(9),3,3)))

# simulate data
n <- 1e5
Z <- rmvnorm(n, sigma = Rho)
U <- pnorm(Z)
X <- matrix(NA, n, 3)
X[,1] <- qf1(U[,1])
X[,2] <- qf2(U[,2])
X[,3] <- qf3(U[,3])

par(mfrow = c(1,3))
invisible(apply(X, 2, hist, breaks = 100, freq = FALSE))

#--- fit the gaussian copula ----------------------------------------------------

gcop <- gcopFit(X = X, n = 512, fitXD = "gc4")

# ok check density estimates
par(mfrow = c(1,3))
for(ii in 1:3) {
  hist(X[,ii], breaks = 100, freq = FALSE)
  curve(dXD(x, xDens = gcop$XDens[[ii]]), add = TRUE, col = "red")
  lines(density(rXD(n = 1e6, gcop$XDens[[ii]])), col = "blue")
}

# check gcop density
# old implementation
source("C:/Users/Jerome/Documents/R/GaussCop/R/old/gausscop-functions.R")

system.time({
  gcop <- gcopFit(X = X, n = 512)
})
system.time({
  gcop2 <- cop.par(X = X, n = 512)
})

# convert gcop2 to gcop format
gcop3 <- lapply(1:3, function(ii) {
  list(ndens = length(gcop2$dens.x[[ii]]), xrng = gcop2$rx[,ii],
       ypdf = gcop2$dens.y[[ii]], ylpdf = gcop2$ldens.y[[ii]],
       ycdf = gcop2$Dens.y[[ii]], mean = gcop2$mean[ii], sd = gcop2$sd[ii])
})

# ok these are not identical because:
# 1. cop.par does normalizes dens.y but not ldens.y
# 2. mean and sd are calculated from sample in cop.par but from density estimate in gcopFit
lapply(1:3, function(ii) {
  gc <- gcop$XDens[[ii]]
  gc2 <- gcop3[[ii]]
  sapply(names(gc), function(jj) range(gc[[jj]] - gc2[[jj]]))
})

# finish conversion
gcop3 <- list(XDens = gcop3, Rho = gcop2$Rho)
class(gcop3) <- "gcop"

dgcop(X = X[1:10,], gCop = gcop3) - dcop(x = X[1:10,], par = gcop2)

#--- speed test for qXD ---------------------------------------------------------

nreps <- 1e6
set.seed(5)
system.time({
  gc1 <- rcop(n = nreps, par = gcop2)
})
set.seed(5)
system.time({
  gc2 <- rgcop(n = nreps, gCop = gcop3)
})

range(gc1 - gc2)

#--- conditional simulation ------------------------------------------------

require(GaussCop)
require(mvtnorm)
# source during debugging phase
source("c:/Users/Jerome/Documents/R/GaussCop/R/rgcopCond.R")

# ok let's start by fitting a multivariate normal
rMn <- function(n, p) {
  if(missing(p)) p <- n
  matrix(rnorm(n*p), n, p)
}

# gives mean and variance
MVCond <- function(Mu, V, xCond, iCond) {
  VOO <- V[iCond,iCond,drop=FALSE]
  VOP <- V[iCond,-iCond,drop=FALSE]
  VPP <- V[-iCond,-iCond,drop=FALSE]
  MuO <- Mu[iCond]
  MuP <- Mu[-iCond]
  # this is not about speed
  A <- t(VOP) %*% solve(VOO)
  list(Mu = MuP + A %*% (xCond - MuO), V = VPP - A %*% VOP)
}

d <- 5
Mu <- rnorm(d)
V <- crossprod(rMn(d))

n <- 1e6
X <- rmvnorm(n, Mu, V)
colnames(X) <- LETTERS[sample(26,d)]
gCop <- gcopFit(X = X, n = 1028)

# univariate
iCond <- sample(d, d-1)
iCond <- rep(TRUE, d)
iCond[sample(d, 1)] <- FALSE
XCond <- X[sample(n,1),iCond]

XMarg <- rgcopCond(n = n, gCop = gCop,
                   XCond = XCond, iCond = iCond, debug = FALSE)

hist(XMarg, breaks = 100, freq = FALSE)
mvc <- MVCond(Mu = Mu, V = V, xCond = XCond, iCond = which(iCond))
curve(dnorm(x, mean = mvc$Mu, sd = sqrt(mvc$V)), col = "red", add = TRUE)

# multivariate
nC <- 3
nM <- d-nC
iCond <- sample(d, nC)
iCond <- rep(TRUE, d)
iCond[sample(d, nM)] <- FALSE
XCond <- X[sample(n,1),iCond]

XMarg <- rgcopCond(n = n, gCop = gCop,
                   XCond = XCond, iCond = iCond)
mvc <- MVCond(Mu = Mu, V = V, xCond = XCond, iCond = which(iCond))

aa <- rnorm(nM)
hist(XMarg %*% aa, breaks = 100, freq = FALSE)
curve(dnorm(x, mean = sum(mvc$Mu * aa),
            sd = sqrt(aa %*% mvc$V %*% aa)), col = "red", add = TRUE)


# check internal multivariate normal functions

require(mvtnorm)

p <- sample(1:10, 1)
n <- 10
sigma <- crossprod(matrix(rnorm(p^2),p))

# doesn't work due to chol pivoting
set.seed(1)
X1 <- rmvnorm(n, sigma = sigma, method = "chol", pre0.9_9994 = TRUE)
set.seed(1)
X2 <- .rmvn(n, sigma = sigma)
range(X1 - X2)

lp1 <- .lmvn(X1, sigma)
lp2 <- dmvnorm(X1, sigma = sigma, log = TRUE)
range(lp1-lp2)

#--- scratch --------------------------------------------------------------------

# ok.  let's do it bit by bit
dcop(x = X[1:10,], par = gcop2, debug = TRUE)
dXD(x = X[1:10,1], xDens = gcop3$XDens[[1]], log = TRUE, debug = TRUE)

# check boundary values
qXD(p = c(-.1, 0, .2, .5, .7, .9, 1, 1.1),
    xDens = gcop3$XDens[[1]], debug = FALSE)
