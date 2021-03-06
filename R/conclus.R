#' Simple wrappers to clustering functions that return just the class memberships
#' @aliases fannyCons hclustCons
#' @usage pamCons(x, k)
#' fannyCons(x, k)
#' hclustCons(x, k)
#' @param x A dissimilarity matrix as produced, for example, by \code{cluster::daisy}
#'   or \code{stats::dist}.
#' @param k The number of clusters to find.
#' @details The \code{conclus} function takes an argument \code{cluster} that should
#'   be a function with precisely 2 arguments, \code{x} and \code{k}. The
#'   function should return an integer vector of cluster memberships for each item
#'   in \code{x}.
#' @export pamCons
pamCons <- function(x, k){
  cluster::pam(x, k, cluster.only=TRUE)
}

#' @export fannyCons
fannyCons <- function(x, k){
  cluster::fanny(x, k, cluster.only=TRUE, diss=TRUE)$clustering
}

#' @export hclustCons
hclustCons <- function(x, k){
  stats::cutree(hclust(x, method="complete"), k)
}

#' Compute connectivity matrix from class memberships
#' @param x Numeric (integer) vector of class memberships.
#' @noRd
connectivity <- function(x){
  k <- length(unique(x)) # number of clusters

  res <- sapply(1:k, function(X) as.numeric(x == X))

  res %*% t(res)
}

#' Compute average connectivity matrix from class memberships
#' @param x The output of \code{subsampleCluster}.
#' @param nms A character vector to be used for row and column names. Defaults
#'   to \code{nms=NULL} and assigning names is not attempted.
#' @noRd
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
#' @noRd
subsampleCluster <- function(diss, k, cluster, subsample, X, seeds=NULL){
  if (!is.null(seeds)){
    set.seed(seeds[[X]])
  }

  res <- rep(0, length.out=nrow(diss))
  names(res) <- rownames(diss)
  if (subsample < 1){ # subsampling
    i <- sample(rownames(diss), size=round(subsample * nrow(diss)))
  } else if (subsample == 1){ # bootstrap sampling
    i <- sample(rownames(diss), size=nrow(diss), replace=TRUE)
  }

  res[i] <- cluster(as.dist(diss[i, i]), k=k)

  cons <- as.numeric(names(res) %in% i)

  list(res, cons %*% t(cons))
}

#' Perform consensus clustering
#' @param diss A dissimilarity matris as returned by, for example, \code{dist}
#'   or \code{daisy}.
#' @param cluster A clustering function that takes 2 arguments: \code{x} and \code{k}
#'   and returns only the class memberships. Functions \code{pamCons} and \code{hclustCons}
#'   are two simple examples. The dissimilarity matrix \code{diss} will be passed
#'   into \code{cluster} so \code{cluster} should NOT coerce \code{x} to be a
#'   dissimilarity matrix.
#' @param subsample The subsampling proportion. Defaults to \code{subsample=0.5}
#'   and 50\% subsampling is performed. If \code{subsample == 1}, bootstrap
#'   sampling (sampling with replacement) is performed.
#' @param K The maximum number of clusters to identify. All values between 2
#'   and \code{K} are used and the consensus clustering matrix returned for
#'   each.
#' @param R The number of random subsamples to run. Defaults to \code{R=100}.
#' @param verbose Whether to report progress. Defaults to \code{verbose=FALSE}.
#' @param ncores The number of cores to use. Defaults to \code{ncores=1} and it
#'   is often the case that \code{conclus} will run faster on a single core than
#'   when it makes the effort to parallelize. To have the function guess the
#'   number of cores, specify \code{ncores=NULL}.
#' @details R random subsamples (or bootstrap samples if \code{subsample = 1}) are
#'   taken from the dissimilarity matrix, and clustering is performed for each value
#'   of k = 2, ..., K. For each value of k, the consensus matrix is computed; each
#'   entry (i, j) represents the average number of times items i and j were in the same
#'   cluster. As such, each element of M is on [0, 1], with 0 or 1 representing perfect
#'   consensus. If items of the concensus matrix are arranged according to cluster
#'   membership, perfect consensus would be represented by a block diagonal form with
#'   blocks full of 1s surrounded by 0s.
#' @return An object of class `conclus'. It contains:
#'   \item{call}{the function call;}
#'   \item{M}{a list, with one element for each k in 2:K, representing the consensus matrices;}
#'   \item{membership}{a matrix with one column, for each k in 2:K, representing the cluster memberships;}
#'   \item{K}{the values of k in 2:K;}
#'   \item{cluster}{the function used to perform the clustering on the subsamples.}
#' @author Harry Southworth
#' @seealso \code{\link{pamCons}},  \code{\link{summary.conclus}}, \code{\link{representatives}}
#' @export conclus
#' @examples
#' # The pluton data
#' cc <- conclus(dist(pluton), K=7) # default PAM clustering
#' ggplot(cc)
#' ggplot(summary(cc))
#' # Do the Gaussian3 and Unform1 examples from Monti et al
#' # First, they used average linkage, so define a new function
#' aveHclustCons <- function(x, k){
#'   stats::cutree(hclust(x, method="average"), k)
#' }
#' # Now pass it into conclus with the Gaussian3 data
#' ccg <- conclus(daisy(Gaussian3), K=6, cluster=aveHclustCons, subsample=.8, R=500, ncores=7)
#' ggplot(ccg, low="white", high="red")
#' s <- summary(ccg)
#' s
#' ggplot(s)
#' # Those are similar to Figures 2 and 3. Do the missing histogram
#' hist(ccg$M[[2]], col="red")
#'
#' # Now Uniform 1
#' ccu <- conclus(daisy(Uniform1), K=6, cluster=aveHclustCons, subsample=.8, R=500, ncores=7)
#' ggplot(ccu, low="white", high="red")
#' su <- summary(ccu)
#' su
#' ggplot(su)
#' hist(c(ccu$M[[2]]), col="green")
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

  # Set up parallelization
  if (is.null(ncores)){
    ncores <- parallel::detectCores() -1
  }
  if (ncores > 1){
    cl <- parallel::makeCluster(ncores)
  }

  for (k in 2:K){
    if (verbose){
      message(paste("Clustering with k =", k))
    }

    if (ncores == 1){
      res <- replicate(R, subsampleCluster(diss, k, cluster, subsample), simplify=FALSE)
    } else {
      seeds <- as.integer(runif(R, -(2^31 - 1), 2^31))

      res <- parallel::parLapply(cl, X=1:R, subsampleCluster,
                                 diss=diss, k=k, cluster=cluster, subsample=subsample)
    }
    M <- averageConnectivity(res, rownames(diss))

    # Reorder M so that items in the same cluster are together
    membership[, k-1] <- cluster(as.dist(1 - M), k=k)
    o <- order(membership[, k-1])

    out[[k - 1]] <- M[o, o]
  } # Close for (k in
  if (ncores > 1){
    parallel::stopCluster(cl)
  }

  colnames(membership) <- paste0("k=", 2:K)
  rownames(membership) <- rownames(diss)

  out <- list(call=theCall, M=out, membership=membership, K=2:K, cluster=cluster)
  class(out) <- "conclus"
  out
}

#' @method ggplot conclus
#' @export
ggplot.conclus <- function(data, mapping, low="peachpuff", high="blue", legend.position="none", ..., environment){
  lvls <- paste0("k=", 2:max(data$K))
  d <- lapply(data$M, as.data.frame)

  d <- lapply(2:(length(d) + 1), function(X) {
    d[[X - 1]] <- mutate(d[[X - 1]], vars = rownames(d[[X - 1]])) %>%
      gather(vars2, M, -vars) %>%
      mutate(k=factor(paste0("k=", X), levels=lvls))
    d[[X - 1]]
  })

  plt <- colorRampPalette(c(low, high))

  p <- list()
  for (i in 1:(max(data$K) - 1)){
    pd <- mutate(d[[i]], vars=factor(vars, levels=unique(vars)),
                 vars2=factor(vars2, levels=unique(vars2)))

    p[[i]] <- ggplot(pd, aes(vars, vars2, fill=M)) +
                geom_tile() +
                scale_fill_gradientn(colors=plt(10)) +
      scale_x_discrete("", breaks=NULL) +
      scale_y_discrete("", breaks=NULL) +
      theme(legend.position=legend.position) +
      facet_wrap(~k)
  }

  do.call("grid.arrange", p)
  invisible(p)
}
