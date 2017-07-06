#' Get cluster and item consensus indices
#' @param x An object of class 'consclus'
#' @param type Character string taking values 'cluster' or 'item'
#' @details If \code{type='cluster'} a matrix with the cluster consensus values
#'   is returned. A \code{print} method exists. If \code{type='item'} an array
#'   is returned with as many rows as items. To access the item consensus for,
#'   say, k=4, use \code{x[,,3]} (remembering that the array extent starts with
#'   k=2). The item consensus values give a measure of belongingness to each item
#'   for each cluster. Prototypes can be extracted by finding the item with the
#'   largest consensus for each cluster.
#' @seealso \code{\link{representatives}}, \code{\link{conclus}}, \code{\link{assignClusters}}
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

#' Get representative items from a consensus clustering
#' @param x The output of a call to \code{conclus}.
#' @param k The number of clusters you believe there to be.
#' @details The function returns the rownames of the representative items, as
#'   judged by them having the highest item consensus within their allocated
#'   cluster.
#' @seealso \code{\link{consensus}}, \code{\link{conclus}}, \code{\link{assignClusters}}
#' @export representatives
representatives <- function(x, k){
  ii <- consensus(x, type="item")
  ii4 <- ii[, 1:k, k-1]

  protos <- apply(ii4, 2, function(x) names(x)[x == max(x)])
  protos
}

#' From conclus object, assign cluster memberships to all items in a dissimilarity matrix
#' @param d A dissimilarity matrix containing a full dataset from which a subset
#'   has been used in \code{conclus}.
#' @param conclus An object of class 'conclus'.
#' @param k The number of clusters you believe there to be.
#' @details With 'large' datasets, dissimilarity matrices get big and it is
#'   sometimes useful to do clustering on the basis of a subset, to conserve
#'   memory. Once clusters have been identified, the cluster members of the items
#'   not in the subset might be required. This function obtains them by finding
#'   the representative items, on the basis of item consensus, and allocating
#'   other items on the basis of their minum distance from the represenetatives.
#' @seealso \code{\link{consensus}}, \code{\link{representatives}}
#' @export assignClusters
assignClusters <- function(d, conclus, k){
  d <- as.matrix(d)
  reps <- representatives(conclus, k=k)

  p_d <- d[, reps]
  apply(p_d, 1, function(X, reps){ reps[X == min(X)] }, reps=names(reps))
}
