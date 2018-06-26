#' Extract a subset of variables from a Gaussian Copula distribution.
#'
#' @param gCop A \code{gaussCop} object (see \code{\link{gcopFit}}).
#' @param subset Integer or logical vector specifying which subset of variables to keep. Default is to keep all variables.
#' @return A \code{gaussCop} object representing the marginal Gaussian Copula distribution on the subset of variables.
#' @examples
#' # simulate data
#' n <- 5e4
#' X <- cbind(rnorm(n, mean = 1, sd = 3),
#'            rnorm(n, mean=4, sd = 0.5),
#'            rt(n, df = 10),
#'            rchisq(n, df = 5),
#'            rnorm(n, mean=10, sd = 10))
#' # fit Gaussian Copula using Kernel method
#' gCop <- gcopFit(X, fitXD = "kernel")
#' # subset gCop
#' isub <- sample(5, 3) # subset
#' gCop.sub <- gcopSub(gCop, subset = isub)
#' # check that subsetted xDensity representations are the same
#' for(ii in 1:3) {
#'   print(identical(gCop$XDens[[isub[ii]]], gCop.sub$XDens[[ii]]))
#' }
#' # check that subsetted correlation matrix is the same
#' print(identical(gCop$Rho[isub,isub], gCop.sub$Rho))
#' @export
gcopSub <- function(gCop, subset) {
  # convert subset to a numeric vector
  nrv <- length(gCop$XDens)
  nm <- names(gCop$XDens)
  if(is.null(nm)) nm <- 1:nrv
  if(missing(subset)) subset <- nm
  iRV <- .getiRV(RVnames = nm, iRV = subset)
  gc <- list(XDens = gCop$XDens[iRV], Rho = gCop$Rho[iRV,iRV,drop=FALSE])
  class(gc) <- "gaussCop"
  gc
}
