#' (Robustly) calculate skewness of a vector.
#' @param x A numeric vector.
#' @param trim A number between 0 and 1 indicating what proportion of values to
#'   trim from the lower and upper ends of \code{x}. Defaults to \code{trim = 0.05}
#'   and 5\% of values are trimmed from each of the lowest and highest values of
#'   \code{x}.
#' @param na.rm What to do about missing values. Defaults to \code{na.rm=TRUE}
#'   and missing values are removed. Otherwise the function will return \code{NA}
#'   if there are missing values in \code{x}.
#' @export
rskew <- function(x, trim=.05, na.rm=TRUE){
  if (trim >= 1 | trim < 0){
    stop("trim should be on [0, 1)")
  }
  if (na.rm){
    x <- x[!is.na(x)]
  }

  x <- trim.hilo(x, trim=trim)

  top <- sum((x - mean(x))^3) / length(x)
  btm <- sd(x)^3

  top / btm
}

#' Trim a proportion of values from the upper and lower extremes of a numeric vector.
#' @param x A numeric vector.
#' @param trim A number on [0, 1).
trim.hilo <- function(x, trim){
  x <- sort(x)
  lo <- floor(length(x) * trim) + 1
  hi <- length(x) + 1 - lo

  lo <- 1:(lo - 1)
  hi <- (hi + 1):length(x)

  x[-c(lo, hi)]
}

#' Robustly scale a numeric vector.
#' @param x A numeric vector.
#' @param na.rm Whether to remove missing values. Defaults to \code{na.rm=TRUE}.
#' @details The vector will be centered by its median and then scaled by its
#'   median absolute deviation. There may be better (more efficient) ways of
#'   doing this (say by using M-estimators of location and scale), but for most
#'   cases, this approach ought to be fine.
#' @export
rscale <- function(x, na.rm=TRUE){
  if (na.rm){
    x <- x[!is.na(x)]
  }

  (x - median(x)) / mad(x)
}

#' (Robustly) deskew a numeric vector.
#' @param x A numeric vector.
#' @param alpha Power parameter in a Box-Cox transformation.
#' @param beta Shift parameter in a Box-Cox transformation.
#' @param na.rm Whether to remove missing values. Defaults to \code{na.rm=TRUE}.
#' @param trim How much to trim from either end of \code{x}. Defaults to \code{trim=0.05}.
#' @return A numeric vector which has been transformed to minimize its (trimmed) skewness.
#'   The vector has attributes \code{par} giving the optimzed values of alpha and beta,
#'   and \code{convergence} which is the convergence code returned by \code{stats::optim}.
#' @details The function trims the data as instructed then uses \code{stats::optim}
#'   to find the Box-Cox parameters that minimize the squared skewness.
#' @export
deskew <- function(x, alpha, beta, na.rm=TRUE, trim=.05){
  if (na.rm){
    x <- x[!is.na(x)]
  }

  x <- trim.hilo(x, trim=trim)

  trans <- function(par=c(a, b), x){
    a <- par[1]
    b <- par[2]
    if (a == 0){
      x <- log(x + b)
    } else {
      x <- ((x + b)^a - 1) / a
    }
  }

  fun <- function(par, x){
    if (any(x + par[2] < 0)){
      10^6
    }
    x <- trans(par, x)
    res <- rskew(x)^2
    if (is.na(res)){
      10^6
    } else {
      res
    }
  }

  res <- optim(par=c(a=0, b= -min(x) + 1), fn=fun, x=x, method="L-BFGS-B",
               lower=c(-Inf, -min(x) + 1))

  out <- trans(res$par, x)
  attr(out, "par") <- res$par
  attr(out, "convergence") <- res$convergence
  out
}
