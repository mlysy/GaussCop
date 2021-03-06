# trimmed mean based not on quantile, but on contribution.
# i.e. eliminate the heaviest 1%.
# if weighted = TRUE, then elimate the top 1% of the weight contributing
# to the sum.
trimmed.mean <- function(x, trim = 0, weighted = FALSE) {
  if(trim == 0) return(mean(x))
  abs.x <- abs(x)
  if(!weighted) q <- quantile(abs.x, probs = 1-trim)
  else q <- wquantile(abs.x, weights = abs.x, probs = 1-trim)
  mean(x[abs.x <= q])
}

# quantile function with weights
wquantile <- function(x, weights = NULL, probs = seq(0, 1, 0.25),
                      names = TRUE) {
  if(is.null(weights)) return(quantile(x, probs = probs, names = names))
  weights <- cumsum(weights[order(x)])
  weights <- weights/max(weights)
  x <- sort(x)
  # for some reason doesn't always work for p = 1
  q <- function(p) {
    if(p == 1) return(max(x))
    x[which.max(weights >= p)]
  }
  ans <- as.numeric(sapply(probs, q))
  if(names) names(ans) <- paste(probs*100, "%", sep = "")
  ans
}

# extract the indices of the conditioning variables
# return the indices as an integer vector, not logical vector.
.getiRV <- function(RVnames, iRV) {
  if(is.character(iRV)) {
    RV <- 1:length(RVnames)
    names(RV) <- RVnames
    iRV <- as.numeric(RV[iRV])
  } else if(is.logical(iRV)) {
    iRV <- which(iRV)
  }
  iRV
}
