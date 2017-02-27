#' Simulation of a Gaussian Copula with Conditioning.
#'
#' @param n Integer number of random draws
#' @param gCop An object of class \code{gcop} describing the joint Gaussian Copula distribution on all variables.
#' @param XCond Vector or matrix; values of the variables on which to condition.
#' @param iCond Vector of integers, logicals, or characters specifying which are the random variables on which to condition.
#' @return An \code{n x (nrv-nCond)} matrix containing the draws from the conditional distribution.
#' @export
rgcopCond <- function(n, gCop, XCond, iCond, debug = FALSE) {
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
  if(debug) browser()
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
  rmvnorm(n, sigma = VPP - IP[(1:nP),,drop=FALSE]) + IP[-(1:nP),,drop=FALSE]
}
