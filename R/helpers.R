#' Robustly scale a numeric vector.
#' @param x A numeric vector.
#' @param na.action What to do when there are missing values. Defaults to
#'   \code{na.action=na.omit} and missing values are returned in the same
#'   places as in \code{x}. The alternative is \code{na.action=na.fail} which
#'   will cause the funtion to fail if there are missing values.
#' @details The vector will be centered by its median and then scaled by its
#'   median absolute deviation. There may be better (more efficient) ways of
#'   doing this (say by using M-estimators of location and scale), but for most
#'   cases, this approach ought to be fine.
#' @export
rscale <- function(x, center=TRUE, scale=TRUE, na.action=na.omit){
  anyNA <- sum(is.na(x))
  if (anyNA){
    n <- length(x)
    x <- na.action(x)
  }

  if (center){
    x <- x - median(x)
  }
  if (scale){
    x <- x / mad(x)
  }

  if (anyNA){
    out <- rep(NA, n)
    out[-c(attributes(x)$na.action)] <- (x - median(x)) / mad(x)
    out
  } else {
    x
  }
}
