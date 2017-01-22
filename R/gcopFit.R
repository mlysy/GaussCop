#' Fit a Gaussian Copula Model.
#'
#' @param X either an \code{N x p} data matrix or a \code{p}-length list of \code{XDens} objects.
#' @param Rho \code{p x p} correlation matrix (i.e., 1's on the diagonal).  Optional if \code{X} is a matrix.
#' @param fitXD string specifying method to fit marginals (see details).
#' @param ... additional arguments to pass to the methods of \code{fitXD} (see details).
#' @return An object of class \code{gcop} (see details).
#' @details A \code{gcop} object is a list with elements:
#' \itemize{
#'   \item \code{XDens}: a list of \code{xDens} objects specifying each marginal distribution (see \code{\link{xDens}}).
#'   \item \code{Rho}: the correlation matrix between the normalized quantiles.
#' }
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
    Z <- sapply(1:ncol(X), function(ii) {
      pXD(q = X[,ii], xDens = XDens[[ii]])
    })
    Rho <- cor(qnorm(Z))
    X <- XDens
  } else if(is.list(X)) {
    if(!all(sapply(X, class) == "xDens")) {
      stop("X must be a matrix or list of xDens objects.")
    }
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
