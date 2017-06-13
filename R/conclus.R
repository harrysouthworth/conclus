#' Simple wrappers to clustering functions that return just the class memberships
#' @aliases fannyCons hclustCons
#' @param x A dissimilarity matrix as produced, for example, by \code{cluster::daisy}
#'   or \code{stats::dist}.
#' @param k The number of clusters to find.
#' @details The \code{conclus} function takes an argument \code{cluster} that should
#'   be a function with precisely to arguments, \code{x} and \code{k}. The
#'   function should return an integer vector of cluster memberships for each item
#'   in \code{x}.
#' @export
pamCons <- function(x, k){
  cluster::pam(x, k, cluster.only=TRUE)
}

#' @export
fannyCons <- function(x, k){
  cluster::fanny(x, k, cluster.only=TRUE, diss=TRUE)$clustering
}

#' @export
hclustCons <- function(x, k){
  stats::cutree(hclust(x, method="complete"), k)
}

#' Compute connectivity matrix from class memberships
#' @param x Numeric (integer) vector of class memberships.
#' @export connectivity
connectivity <- function(x){
  k <- length(unique(x)) # number of clusters

  res <- sapply(1:k, function(X) as.numeric(x == X))

  res %*% t(res)
}

#' Compute average connectivity matrix from class memberships
#' @param x The output of \code{subsampleCluster}.
#' @param nms A character vector to be used for row and column names. Defaults
#'   to \code{nms=NULL} and assigning names is not attempted.
#' @export averageConnectivity
averageConnectivity <- function(x, nms=NULL){
  # Compute the average connectivity matrix
  M <- apply(sapply(x, function(X) connectivity(X[[1]]), simplify="array"), 1:2, sum)

  # Get number of times each item was included in a subsample
  N <- apply(sapply(x, function(X) X[[2]], simplify="array"), 1:2, sum)

  M <- M / N

  if (any(N == 0)){ # deal wtih 0 / 0
    M[N == 0] <- 0
  }

  checkM(M)

  # Set dim names and rescale
  if (!is.null(nms)){
    colnames(M) <- rownames(M) <- nms
  }

  M
}

#' Perform clustering on a subsample of a dissimilarity matrix
#' @param diss A dissimilarity matrix.
#' @param k The number of clusters to find.
#' @param cluster A function that returns a character vector of cluster memberships.
#' @param subsample A number on (0, 1) dictating the proportion of items to subsample.
#' @param X Arugment used by calls to \code{parLapply}, not to be specified by the user.
#' @param seeds Vector of seeds used by parLapply
#' @details This function is called by conclus and is split out separately
#'   purely to enable unit testing
subsampleCluster <- function(diss, k, cluster, subsample, X, seeds=NULL){
  if (!is.null(seeds)){
    set.seed(seeds[[X]])
  }

  res <- rep(0, length.out=nrow(diss))
  names(res) <- rownames(diss)
  i <- sample(rownames(diss), size=round(subsample * nrow(diss)))

  res[i] <- cluster(as.dist(diss[i, i]), k=k)

  cons <- as.numeric(names(res) %in% i)

  list(res, cons %*% t(cons))
}

#' Perform consensus clustering
#' @param diss A dissimilarity matris as returned by, for example, \code{dist}
#'   or \code{daisy}. It will be immediately turned into a matrix via \code{as.matrix}.
#' @param cluster A clustering function that takes 2 arguments: \code{x} and \code{k}
#'   and returns only the class memberships. Functions \code{pamCons} and \code{hclustCons}
#'   are two simple examples. The dissimilarity matrix \code{diss} will be passed
#'   into \code{cluster} so \code{cluster} should NOT coerce \code{x} to be a
#'   dissimilarity matrix.
#' @param subsample The subsampling proportion. Defaults to \code{subsample=0.5}
#'   and 50\% subsampling is performed.
#' @param K The maximum number of clusters to identify. All values between 2
#'   and \code{K} are used and the consensus clustering matrix returned for
#'   each.
#' @param R The number of random subsamples to run. Defaults to \code{R=100}.
#' @param verbose Whether to report progress. Defaults to \code{verbose=FALSE}.
#' @param ncores The number of cores to use. Defaults to \code{ncores=1} and it
#'   is often the case that \code{conclus} will run faster on a single core than
#'   when it makes the effort to parallelize. To have the function guess the
#'   number of cores, specify \code{ncores=NULL}.
#' @export conclus
conclus <- function(diss, cluster=pamCons, subsample=.5, K=NULL, R=100, verbose=FALSE,
                    ncores=1){
  theCall <- match.call()

  checkDiss(diss)
  checkK(K)
  checkSubsample(subsample)
  checkClusterFun(cluster)

  diss <- as.matrix(diss) # Need to to this for subsampling

  out <- list()
  membership <- matrix(0, nrow=nrow(diss), ncol=K-1)
  rownames(membership) <- rownames(diss)

  for (k in 2:K){
    if (verbose){
      message(paste("Clustering with k =", k))
    }

    # Set up parallelization
    if (is.null(ncores)){
      ncores <- parallel::detectCores() -1
    }
    if (ncores == 1){
      res <- replicate(R, subsampleCluster(diss, k, cluster, subsample), simplify=FALSE)
    } else {
      cl <- parallel::makeCluster(ncores)
      seeds <- as.integer(runif(R, -(2^31 - 1), 2^31))

      res <- parallel::parLapply(cl, X=1:R, subsampleCluster,
                                 diss=diss, k=k, cluster=cluster, subsample=subsample)
      parallel::stopCluster(cl)
    }
    M <- averageConnectivity(res, rownames(diss))

    # Reorder M so that items in the same cluster are together
    membership[, k-1] <- cluster(as.dist(1 - M), k=k)
    o <- order(membership[, k-1])

    out[[k - 1]] <- M[o, o]
  }

  colnames(membership) <- paste0("k=", 2:K)
  rownames(membership) <- rownames(diss)

  out <- list(M=out, membership=membership, K=2:K, cluster=cluster, call=theCall)
  class(out) <- "conclus"
  out
}

#' @export
ggplot.conclus <- function(data, mapping, low="white", high="blue", legend.position="none", ..., environment){
  lvls <- paste0("k=", 2:max(data$K))
  d <- lapply(data$M, as.data.frame)

  d <- lapply(2:(length(d) + 1), function(X) {
    d[[X - 1]] <- mutate(d[[X - 1]], vars = rownames(d[[X - 1]])) %>%
      gather(vars2, M, -vars) %>%
      mutate(k=factor(paste0("k=", X), levels=lvls))
    d[[X - 1]]
  })


  p <- list()
  for (i in 1:(max(data$K) - 1)){
    pd <- mutate(d[[i]], vars=factor(vars, levels=unique(vars)),
                 vars2=factor(vars2, levels=unique(vars2)))
    p[[i]] <- ggplot(pd, aes(vars, vars2)) +
                geom_tile(aes(fill=M), color=low) +
                scale_fill_gradient(low=low, high=high) +
      scale_x_discrete("", breaks=NULL) +
      scale_y_discrete("", breaks=NULL) +
      theme(legend.position=legend.position) +
      facet_wrap(~k)
  }

  do.call("grid.arrange", p)
}
