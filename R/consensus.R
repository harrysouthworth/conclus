#' Get cluster and item consensus indices
#' @param x An object of class 'consclus'
#' @param type Character string taking values 'cluster' or 'item'
#' @export consensus
consensus <- function(x, type="cluster"){
  if (class(x) != "conclus"){
    stop("x should have class 'conclus'")
  }

  if (type == "cluster"){
    clusterConsensus(x)
  } else if (type == "item"){
    itemConsensus(x)
  } else {
    stop("type should be either 'cluster' or 'item'")
  }
}

clusterConsensus <- function(x){
  cons <- matrix(0, nrow=max(x$K - 1), ncol=max(x$K))
  for (k in x$K){
    o <- apply(x$membership, 2, order)
    for (i in 1:k){
      # Each item of x$M is ordered, x$membership is not
      oi <- x$membership[o[, k-1], k-1] == i
      subM <- x$M[[k-1]][oi, oi, drop=FALSE]
      diag(subM) <- 0

      N <- nrow(subM)
      if (is.null(N)){
        # A single item in a cluster
        message(paste0("With k=", k, " cluster ", i, " has a single member."))
        cons[k-1, i] <- 1
      } else { # Mean of whole thing is same as mean of lower triangle
        cons[k-1, i] <- mean(subM)
      }
    }
  }

  rownames(cons) <- paste0("k=", x$K)
  colnames(cons) <- paste("cluster", 1:ncol(cons))
  t(cons)
}

itemConsensus <- function(x){
  cons <- array(0, dim = c(nrow(x$M[[1]]), max(x$K), max(x$K -1)))

  for (k in x$K){
    M <- x$M[[k-1]]
    diag(M) <- 0

    membership <- x$membership
    o <- order(membership[, k-1])
    membership <- membership[o, k-1]

    Nk <- table(membership)
    cN <- c(0, cumsum(Nk))

    for (i in 1:k){
      for (j in 1:nrow(M)){
        if (membership[j] == i){
          if (Nk[i] == 1){
            message(paste0("With k=", k, " cluster ", i, " has a single member."))
            cons[j, i, k-1] <- 1
          } else {
            cons[j, i, k-1] <- sum(M[j, (cN[i] + 1):cN[i+1]], diag=FALSE) / (Nk[i] - 1)
          }
        } else {
          cons[j, i, k-1] <- sum(M[j, (cN[i] + 1):cN[i+1]]) / Nk[i]
        }
      }

      # Need to undo the result of using order(x$membership[, k-1]) above
      cons[o, i, k-1] <- cons[, i, k-1]
    }
  }

  dimnames(cons) <-  list(rownames(x$membership),
                          paste("cluster", 1:max(x$K)),
                          paste0("k=", x$K))
  cons
}
