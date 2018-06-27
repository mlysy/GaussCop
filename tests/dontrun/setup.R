#--- setup file for GaussCop package --------------------------------------------

require(devtools)

pkg.path <- "C:/Users/Jerome/Documents/R/GaussCop"
build.path <- "C:/Users/Jerome/Documents/R/build"

document(pkg = pkg.path)
install(pkg = pkg.path)

build(pkg = pkg.path, path = build.path)

#--- package structure ----------------------------------------------------------

# the package does the following things:
# 1.  extended density estimation.  so an xdensity object contains a 2-column matrix of density evaluations on a grid, and something to put outside of the grid.
# for now this is some normal which tries to smoothly extend the grid-based estimator.
# 2.  the grid can be calculated by kernel smoothing, or using a Box-Cox + Gram-Charlier approximation, which works well when the underlying density isn't too far from normal.
# 3.  then provides simple algorithms for gaussian copula density and random sampling.

# ok: function for density estimation
# very simple, but drawback is that we don't see arguments to gc4.fit
xdensity <- function(x, n = 512, from, to, mean, sd, any0 = FALSE,
                     method = c("kernel", "gc4"), ...) {}

# ok now arguments to these are far more clear.
kernelXD <- function(x, n = 512, from, to, mean, sd, any0 = FALSE, ...) {}
gc4XD <- function(x, n = 512, from, to, mean, sd, any0 = FALSE,
                  lambda = NA, alpha = 0, trim = .01, ...) {}

# now this extends a density matrix representation to an xDens object.
xdensity <- function(XY, mean, sd) {}

# these will return objects of type xDens which contains:
# xrng, xd, (dens.x will be depreciated), ypdf, ylpdf, ycdf

# ok define these
dXD <- function(x, xDens, log = FALSE) {}
qXD <- function(p, xDens, log = FALSE) {}
pXD <- function(q, xDens, log = FALSE) {}
rXD <- function(n, xDens, log = FALSE) {}

# to calculate the correlation matrix
# XDens is a list of xDens objects.
RhoXD <- function(X, XDens) {}

gcopFit <- function(X, Rho, methodXD = c("kernel", "gc4"), ...) {}

# so have underlying call to gc4.fit
