#' @title Fit a Gaussian Copula model.
#' @description Fit a Gaussian Copula model to be used for simulation and density calculation. 
#' @param X Either an \code{N x p} data matrix or a \code{p}-length list of \code{XDens} objects.
#' @param Rho \code{p x p} correlation matrix (i.e., 1's on the diagonal).  Optional if \code{X} is a matrix.
#' @param fitXD String specifying method to fit marginals (see \code{\link{XDens}}).
#' @param ... Additional arguments to pass to the methods of \code{fitXD} (see \code{\link{XDens}}).
#' @return An object of class \code{gcop}, i.e., a list with elements:
#' \describe{
#'   \item{\code{XDens}}{A list of \code{xDens} objects specifying each marginal distribution (see \code{\link{xDens}}).}
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
gcopFit <- function(X, Rho, fitXD = c("kernel", "gc4", "normal"), ...) {
  if(is.matrix(X)) {
    # marginal densities
    fitXD <- match.arg(fitXD)
    if(fitXD == "kernel") {
      XDens <- apply(X, 2, kernelXD, ...)
    } else if(fitXD == "gc4") {
      XDens <- apply(X, 2, gc4XD, ...)
    } else if(fitXD == "normal") {
      XDens <- apply(X, 2, normalXD, ...)
    }
    # correlation
    U <- sapply(1:ncol(X), function(ii) {
      pXD(q = X[,ii], xDens = XDens[[ii]])
    })
    if(missing(Rho)) Rho <- cor(qnorm(U))
    X <- XDens
  } else if(is.list(X)) {
    if(!all(sapply(X, class) == "xDens")) {
      stop("X must be a matrix or list of xDens objects.")
    }
    if(missing(Rho)) stop("X is not a matrix so Rho must be supplied.")
  } else {
    stop("X must be a matrix or list of xDens objects.")
  }
  # create gcop object
  nrv <- length(X)
  rvnames <- names(X)
  colnames(Rho) <- rvnames
  rownames(Rho) <- rvnames
  gcp <- list(XDens = X, Rho = Rho)
  class(gcp) <- "gcop"
  gcp
}
