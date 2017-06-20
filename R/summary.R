# Get cumulative distribution functions and areas under them for a consclus object
cdf.conclus <- function(x){
  btm <- nrow(x$M[[1]]) * (nrow(x$M[[1]]) - 1) / 2

  res <- list()

  ii <- lower.tri(x$M[[1]], diag=FALSE)

  for (i in 1:length(x$M)){
    M <- c(x$M[[i]][ii])

    cc <- sort(unique(c(0, M)))

    cdf <- sapply(cc, function(X) sum(M <= X)) / btm

    cdf <- as.data.frame(cbind(index=cc, cdf=cdf))
    cdf$k <- paste0(x$K[i])

    res[[i]] <- cdf
  }

  res <- bind_rows(res)

  names(res) <- c("index", "CDF", "k")

  res
}


#' Compute and return summaries of a 'conclus' object.
#' @param object An object of class 'conclus'.
#' @param ... Not used.
#' @details The function computes the distribution function of the elements of the
#'   consensus matrices returned by the call to \code{conclus}. From these, it
#'   then returns the ares under the curves, and the incremental proportional changes to those areas as
#'   the number of clusters k is increeased. It also computes the 'cluster consensus'
#'   which enables the user to see which of the clusters identified are the most
#'   stable. These statistics are intended to aid the user in selecting how many
#'   clusters they might reasonably believe to exist. They should be used in
#'   conjunction with the plot of the consensus matrices, generated via
#'   \code{ggplot(objec t)}, for the object that is passed to the \code{summary}
#'   function.
#' @export
summary.conclus <- function(object, ...){
  cdf <- cdf.conclus(object)
  auc <- auc.conclus(cdf)
  delta <- delta.conclus(auc)
  cons <- consensus(object, type="cluster")

  res <- list(call=object$call, CDF=cdf, AUC=auc, Delta=delta, Consensus=cons)

  class(res) <- "summary.conclus"
  res
}

auc.conclus <- function(x){
  auc <- vector(length=length(unique(x$k)))

  for (i in 1:length(auc)){
    d <- x[x$k == i + 1, ]
    xx <- d$index
    yy <- d$CDF

    index <- 2:length(xx)

    auc[i] <- (xx[index] - xx[index - 1]) %*% (yy[index] + yy[index - 1]) / 2
  }

  auc
}

delta.conclus <- function(auc){
  # auc should be the output of auc.conclus
  # Need to account for cases where AUC can decline as number of clusters increases
  # See the textual description in Monti et al, beneath the formulae that are
  #   given for delta
  delta <- rep(0, length(auc))
  mauc <- delta[1] <- auc[1]
  for (i in 2:length(auc)){
    delta[i] <- (auc[i] - mauc) / mauc
    mauc <- max(auc[i], mauc)
  }
  delta
}

#' @method print conclus
#' @export
print.conclus <- function(x, digits=3, ...){
  print(x$call)
  x <- summary(x)
  cat("\nAUC\n")
  cat(round(x$AUC, digits=digits))
  cat("\nDelta\n")
  cat(round(x$Delta, digits=digits))
  invisible()
}

#' @method print summary.conclus
#' @export
print.summary.conclus <- function(x, digits=3, ...){
  print(x$call)

  dd <- matrix(round(c(x$AUC, x$Delta), digits=digits), nrow=2, byrow=TRUE)
  rownames(dd) <- c("AUC", "Delta")
  colnames(dd) <- colnames(x$Consensus)
  cat("\nAUC & Delta\n")
  print(dd)
  cat("\nCluster consensus\n")
  con <- format(x$Consensus, zero.print=FALSE, digits=digits)
  print(noquote(con))
  invisible()
}

#' @method ggplot summary.conclus
#' @export
ggplot.summary.conclus <- function(data, mapping=NULL, legend.position="bottom", ..., environment){
  auc <- data$AUC
  delta <- data$Delta
  data <- data$CDF

  lvls <- paste0("k=", unique(data$k))

  data$k <- factor(data$k, labels=lvls)

  cdfplot <-
  ggplot(data, aes(index, CDF, color=k)) +
    geom_step(size=2, alpha=.7) +
    theme(legend.title=element_blank(), legend.position=legend.position) +
    scale_x_continuous("Consensus index", limits=c(0, 1)) +
    scale_y_continuous("CDF", limits=c(0, 1))

  k <- rep(2:(length(lvls) + 1), 2)
  value <- c(auc, delta)
  what <- rep(c("AUC", "Delta(k)"), each=length(lvls))

  auc <- data.frame(k=k, value=value, what=what)

  aucplot <-
  ggplot(auc, aes(k, value)) +
    geom_point(size=3, color="blue") +
    geom_line(color="blue", size=1.5) +
    scale_x_continuous("Number of clusters") +
    scale_y_continuous("Value") +
    facet_wrap(~ what, scales="free_y", labeller=label_parsed, ncol=1)

  grid.arrange(cdfplot, aucplot, ncol=2)

  invisible(list(cdfplot, aucplot))
}
