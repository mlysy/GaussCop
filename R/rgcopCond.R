#' Conditional simulation from a Gaussian Copula distribution.
#'
#' @param n Integer number of random draws.
#' @param gCop An object of class \code{gcop} describing the joint Gaussian Copula distribution on all variables.
#' @param XCond Vector or matrix; values of the variables on which to condition.
#' @param iCond Vector of logicals, or characters specifying which are the random variables on which to condition.
#' @return An \code{n x (nrv-nCond)} matrix containing the draws from the conditional distribution.
#' @examples 
#' require(GaussCop)
#' #conditionally simulate data and plot it
#' n = 1e4
#' dat1 = rnorm(n, mean = 1, sd = 3)
#' dat = cbind(dat1,
#'             rnorm(n, mean = dat1, sd = 2),
#'             rnorm(n, mean=dat1, sd = 10))
#' 
#' plot(dat[, c(2, 3)]) # plot a subset of the data
#' 
#' # fit Gaussian Copula using gc4 method
#' temp.cop = gcopFit(X = dat, fitXD = "gc4")
#' 
#' # simulate data from Copula model and add it to plot, should blend in
#' new.data = rgcopCond(n, gCop = temp.cop, XCond = as.matrix(dat1), iCond = !c(FALSE, TRUE, TRUE))
#' points(new.data, cex = 0.5, col="red") # plot points from subsetted copula
#' @export
rgcopCond <- function(n, gCop, XCond, iCond) {
  if(class(gCop) != "gcop") stop("gCop must be a gcop object.")
  # determine conditioning variables
  nrv <- length(gCop$XDens)
  nm <- names(gCop$XDens)
  if(!is.matrix(XCond)) XCond <- t(XCond)
  if(missing(iCond)) iCond <- colnames(XCond)
  iCond <- .getiRV(RVnames = nm, iRV = iCond)
  nCond <- length(iCond)
  ## iCond <- .getiCond(iCond, XCond, nm)
  ## nCond <- length(iCond)
  ## if(!is.matrix(XCond)) {
  ##   XCond <- matrix(XCond, ncol = nCond)
  ## }
  
  if(!all(iCond %in% 1:nrv) || nCond != ncol(XCond)) {
    stop("Incorrect specification of iCond.")
  }
  # map conditioning variables to desired normals
  for(ii in 1:nCond) {
    jj <- iCond[ii]
    XCond[,ii] <- qnorm(pXD(q = XCond[,ii], xDens = gCop$XDens[[jj]]))
  }
  # generate conditional normals
  iMarg <- (1:nrv)[-iCond]
  XMarg <- .rmvnCond(n, XO = XCond,
                     VOO = gCop$Rho[iCond,iCond,drop=FALSE],
                     VPP = gCop$Rho[iMarg,iMarg,drop=FALSE],
                     VOP = gCop$Rho[iCond,iMarg,drop=FALSE])
  colnames(XMarg) <- nm[iMarg]
  # invert cdfs to obtain marginals
  for(ii in 1:(nrv-nCond)) {
    jj <- iMarg[ii]
    XMarg[,ii] <- qXD(p = pnorm(XMarg[,ii]), xDens = gCop$XDens[[jj]])
  }
  XMarg
}

# extract the indices of the conditioning variables
# return the indices as an integer vector, not logical vector
.getiCond <- function(iCond, XCond, rv.names) {
  nrv <- length(rv.names)
  if(missing(iCond)) {
    if(is.matrix(XCond)) {
      iCond <- colnames(XCond)
    } else {
      iCond <- names(XCond)
    }
  }
  if(is.character(iCond)) {
    x <- 1:nrv
    names(x) <- rv.names
    iCond <- as.numeric(x[iCond])
  } else if(is.logical(iCond)) {
    iCond <- which(iCond)
  }
  iCond
}

.rmvnCond <- function(n, XO, VOO, VPP, VOP) {
  nP <- nrow(VPP)
  IP <- rbind(t(VOP), XO) %*% solve(VOO, VOP)
  .rmvn(n, sigma = VPP - IP[(1:nP),,drop=FALSE]) + IP[-(1:nP),,drop=FALSE]
}
