#' Marginal subset of a gcop object.
#'
#' @param gCop A \code{gcop} object.
#' @param subset Vector specifying which subset of variables to keep.  Default is all.
#' @return A \code{gcop} with variables possibly removed.
#' @export
gcopSub <- function(gCop, subset) {
  # convert subset to a numeric vector
  nrv <- length(gCop$XDens)
  nm <- names(gCop$XDens)
  if(is.null(nm)) nm <- 1:nrv
  if(missing(subset)) subset <- nm
  iRV <- .getiRV(RVnames = nm, iRV = subset)
  gc <- list(XDens = gCop$XDens[iRV], Rho = gCop$Rho[iRV,iRV,drop=FALSE])
  class(gc) <- "gcop"
  gc
}
