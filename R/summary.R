# Get cumulative distribution functions and areas under them for a consclus object
cdf.conclus <- function(x){
  btm <- nrow(x$M[[1]]) * (nrow(x$M[[1]]) - 1) / 2

  res <- list()

  ii <- lower.tri(x$M[[1]], diag=FALSE)

  for (i in 1:length(x$M)){
    M <- c(x$M[[i]][ii])
    #M <- M[M > 0]
    cc <- sort(M)

    cdf <- sapply(cc, function(X) sum(M <= X)) / btm

    cdf <- as.data.frame(cbind(index=cc, cdf=cdf))
    cdf$k <- paste0(x$K[i])

    res[[i]] <- cdf
  }

  res <- bind_rows(res)

  names(res) <- c("index", "CDF", "k")

  res
}

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
    w <- diff(d$index)
    auc[i] <- sum(w * d$CDF[-1])
  }

  auc
}

delta.conclus <- function(auc){
  # auc should be the output of auc.conclus
  # Need to account for cases where AUC can decline as number of clusters increases
  delta <- mauc <- auc[1]
  for (i in 2:length(auc)){
    delta <- c(delta, (auc[i] - mauc) / mauc)
    mauc <- max(auc[i], mauc)
  }
  delta
}


#' @export
plot.summary.conclus <- function(x, ...){
  auc <- x$AUC
  x <- x$CDF
  plot(x[, 1:2], type="l", lwd=2, xlab="Consensus index", ylab="CDF")
  for (i in 3:ncol(x)){
    lines(x[, c(1, i)], col=i+1, lwd=2)
  }

  invisible()
}

#' @export
print.conclus <- function(x, digits=3, ...){
  print(x$call)
  x <- summary(x)
  cat("\nAUC\n")
  print(x$AUC, digits=digits)
  cat("\nDelta\n")
  print(x$Delta, digits=digits)
  invisible()
}

#' @export
print.summary.conclus <- function(x, digits=3, ...){
  print(x$call)
  cat("\nAUC\n")
  print(x$AUC, digits=digits)
  cat("\nDelta\n")
  print(x$Delta, digits=digits)
  cat("\nCluster consensus\n")
  con <- format(x$Consensus, zero.print=FALSE, digits=digits)
  print(noquote(con))
  invisible()
}

#' @export
ggplot.summary.conclus <- function(data, mapping=NULL, legend.position="bottom", ..., environment){
  delta <- data$Delta
  data <- data$CDF

  lvls <- paste0("k=", unique(data$k))

  data$k <- factor(data$k, labels=lvls)

  cdfplot <-
  ggplot(data, aes(index, CDF, color=k)) +
    geom_line(size=2) +
    theme(legend.title=element_blank(), legend.position=legend.position) +
    scale_x_continuous("Consensus index", limits=c(0, 1)) +
    scale_y_continuous("CDF", limits=c(0, 1))

  delta <- data.frame(k=2:(length(lvls) + 1), delta=delta)

  aucplot <-
  ggplot(delta, aes(k, delta)) +
    geom_point(size=3, color="blue") +
    geom_line(color="blue", size=1.5) +
    scale_x_continuous("Number of clusters") +
    scale_y_continuous(expression(Delta(AUC)))

  grid.arrange(cdfplot, aucplot, ncol=2)

  invisible(list(cdfplot, aucplot))
}
