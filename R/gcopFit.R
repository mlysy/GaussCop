#' Fit a Gaussian Copula model.
#'
#' @param X Either an \code{N x p} data matrix or a \code{p}-length list of \code{XDens} objects.
#' @param Rho \code{p x p} correlation matrix (i.e., 1's on the diagonal).  Optional if \code{X} is a matrix.
#' @param fitXD String specifying method to fit marginals (see \code{\link{xDensity}}).
#' @param ... Additional arguments to pass to the methods of \code{fitXD} (see \code{\link{xDensity}}).
#' @return An object of class \code{gaussCop}, i.e., a list with elements:
#' \describe{
#'   \item{\code{XDens}}{A list of \code{xDensity} objects specifying each marginal distribution (see \code{\link{xDensity}}).}
#'   \item{\code{Rho}}{The correlation matrix between the normalized quantiles.}
#' }
#' @examples
#' # simulate data and plot it
#' n = 5e4
#' dat = cbind(rnorm(n, mean = 1, sd = 3),
#'             rnorm(n, mean=4, sd = 0.5))
#' plot(dat, cex=0.5)
#' # fit Gaussian Copula using Kernel method
#' temp.cop = gcopFit(X = dat, fitXD = "kernel")
#' # simulate data from Copula model and add it to plot, should blend in
#' new.data = rgcop(100, temp.cop)
#' points(new.data, cex = 0.5, col="red")
#' @export
gcopFit <- function(X, Rho, fitXD = c("kernel", "gc4"), ...) {
  if(is.matrix(X)) {
    # marginal densities
    fitXD <- match.arg(fitXD)
    if(fitXD == "kernel") {
      XDens <- apply(X, 2, kernelXD, ...)
    } else if(fitXD == "gc4") {
      XDens <- apply(X, 2, gc4XD, ...)
    }
    # depreciated, as can be set by matrixXD
    ## else if(fitXD == "normal") {
    ##   XDens <- apply(X, 2, normalXD, ...)
    ## }
    # correlation
    U <- sapply(1:ncol(X), function(ii) {
      pXD(q = X[,ii], xDens = XDens[[ii]])
    })
    if(missing(Rho)) Rho <- cor(qnorm(U))
    X <- XDens
  } else if(is.list(X)) {
    if(!all(sapply(X, class) == "xDensity")) {
      stop("X must be a matrix or list of xDensity objects.")
    }
    if(missing(Rho)) stop("X is not a matrix so Rho must be supplied.")
  } else {
    stop("X must be a matrix or list of xDensity objects.")
  }
  # create gcop object
  nrv <- length(X)
  rvnames <- names(X)
  colnames(Rho) <- rvnames
  rownames(Rho) <- rvnames
  gcp <- list(XDens = X, Rho = Rho)
  class(gcp) <- "gaussCop"
  gcp
}
