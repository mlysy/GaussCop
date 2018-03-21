#' Extract a subset of variables from a Gaussian Copula distribution.
#'
#' @param gCop A \code{gaussCop} object (see \code{\link{gcopFit}}).
#' @param subset Integer or logical vector specifying which subset of variables to keep. Default is to keep all variables.
#' @return A \code{gaussCop} object representing the marginal Gaussian Copula distribution on the subset of variables.
#' @examples
#' # simulate data and plot it
#' n <- 5e4
#' dat <- cbind(rnorm(n, mean = 1, sd = 3),
#'              rnorm(n, mean=4, sd = 0.5),
#'              rt(n, df = 5),
#'              rchisq(n, df = 3),
#'              rnorm(n, mean=10, sd = 10))
#' plot(dat[, c(3, 4)]) # plot a subset of the data
#' # fit Gaussian Copula using Kernel method
#' temp.cop <- gcopFit(X = dat, fitXD = "gc4")
#' sub.cop <- gcopSub(temp.cop, c(3, 4)) # subset copula
#' # simulate data from Copula model and add it to plot, should blend in
#' new.data <- rgcop(100, sub.cop)
#' points(new.data, cex = 0.5, col="red") # plot points from subsetted copula
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
