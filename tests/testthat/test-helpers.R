context("Helper functions")

test_that("checkM fails when it should", {
  m <- matrix(1:9, ncol=3)
  expect_error(checkM(m), label="checkM: asymmetric matrix causes error")

  m <- diag(rep(1, 3))
  m[1, 2] <- m[2, 1] <- NA
  expect_error(checkM(m), label="checkM: NAs in matrix not allowed")

  m <- diag(1:3)
  expect_error(checkM(m), label="checkM: maximum value of matrix not > 1")

  m <- diag(1 / 2:4)
  expect_error(checkM(m), label="checkM: maximum value of matarix not < 1")

  m <- diag(rep(1, 3))
  m[1, 2] <- m[2, 1] <- -.1
  expect_error(checkM(m), label="checkM: fails when there are negatives")

  m <- diag(c(1, .5, 1))
  expect_error(checkM(m), label="checkM: fails when diagonals not all 1")
})

test_that("cluster functions return correct number of clusters", {
  for (k in 2:5){
    x <- pamCons(daisy(pluton), k=k)
    expect_equal(min(x), 1, label="pamCons: minimum cluster label is 1")
    expect_equal(max(x), k, label=paste("pamCons: maixmum cluster label is", k))

    x <- hclustCons(daisy(pluton), k=k)
    expect_equal(min(x), 1, label="hclustCons: minimum cluster label is 1")
    expect_equal(max(x), k, label=paste("hclustCons: maixmum cluster label is", k))

    x <- fannyCons(daisy(pluton), k=k)
    expect_equal(min(x), 1, label="fannyCons: minimum cluster label is 1")
    expect_equal(max(x), k, label=paste("fannyCons: maixmum cluster label is", k))
  }
})

test_that("inappropriate cluster function causes failure", {
  cf <- function(x, k){
    fanny(x, k=k, diss=TRUE)
  }
  expect_error(conclus(daisy(pluton), K=5, cluster=cf), label="checkClusterFun: returning non-integer vector causes error")

  cf <- function(x, k, z=NULL){
    pam(x, k=k, cluster.only=TRUE)
  }
  expect_error(conclus(daisy(pluton), K=5, cluster=cf), label="checkClusterFun: function with wrong arguments causes error")
})

test_that("connectivity behaves as expected", {
  for (i in 1:10){
    x <- sample(1:5, size=100, replace=TRUE)

    co <- connectivity(x)

    expect_equal(nrow(co), ncol(co), label="connectivity: result is square matrix")
    expect_equal(nrow(co), length(x), label="connectivity: result has correct dimensions")
    expect_true(all(diag(co) == 1), label="connectivity: diagonals are all 1")
    expect_true(isSymmetric(co), label="connectivity: output is symmetric")
    expect_true(all(sort(unique(c(co))) == c(0:1)), label="connectivity: all entries are 0 or 1")

    # All columns in same class should be identical (same for rows, symmetry tested above)
    for (j in unique(x)){
      class <- co[, x == j]
      expect_true(all(apply(class, 2, function(X) all(X == class[, 1]))),
                  label="connectivity: all same class membership rows and columns are identical")
    }
  }
})

test_that("clusters on subsamples are as expected", {
  for (i in 1:10){
    k <- sample(2:5, size=1)
    subsample <- runif(1, .5, .9)

    nsamp <- round(nrow(pluton) * (1 - subsample))

    ss <- subsampleCluster(as.matrix(daisy(pluton)), k=k, cluster=pamCons, subsample=subsample)

    # ss[[1]] should be a vector of class memberships, 0 where item was not in sample
    expect_equal(min(ss[[1]]), 0, label="subsampleCluster: minimum class membership is 0")
    expect_equal(max(ss[[1]]), k, label="paste(subsampleCluster: maximum class membership is", k)
    expect_equal(sum(ss[[1]] == 0), nsamp,
                 label=paste("subsampleCluster: number of 0s is nrow(data) * round(1 -", subsample, ")"))

    # ss[[2]] should be binary matrix indicating if items were in sample
    expect_true(all(sort(unique(c(ss[[2]]))) == 0:1),
                label="subsampleCluster: inclusion matrix is binary")
    expect_equal(unique(dim(ss[[2]])), nrow(pluton),
                 label="subsampleCluster: inclusion matrix has correct dimensions")
    expect_true(isSymmetric(ss[[2]]), label="subsampleCluster: inclusion matrix is symmetric")

    ucs <- unique(colSums(ss[[2]]))
    expect_equal(sort(ucs), c(0, nrow(pluton) - nsamp),
                 label="subsampleCluster: inclusion matrix has correct numbers of 0s and 1s")

    # Get indices of excluded items off ss[[1]]
    excluded <- (1:nrow(pluton))[ss[[1]] == 0]

    expect_equal(apply(ss[[2]][, excluded], 2, sum), rep(0, nsamp),
                 label="subsampleCluster: excluded items from memberships match 0s in inclusion matrix")
  } # close for (i in 1:)
})

test_that("averaging connectivity matrices behaves as expected", {
  # We've already tested connectivity(.) so this is largely a test of the calls
  # to apply and the construction of N within averageConnectivity

  for (i in 1:10){
    k <- sample(2:5, size=1)
    subsample <- runif(1, .5, .9)
    R <- sample(5:20, size=1)

    x <- replicate(R, subsampleCluster(as.matrix(daisy(pluton)), k=k, cluster=pamCons, subsample=subsample),
                   simplify=FALSE)

    totConn <- totN <- matrix(0, nrow=nrow(pluton), ncol=nrow(pluton))
    rownames(totConn) <- colnames(totConn) <- rownames(totN) <- colnames(totN) <- rownames(pluton)

    for (j in 1:R){
      totConn <- totConn + connectivity(x[[j]][[1]])
      totN <- totN + x[[j]][[2]]
    }

    ac <- averageConnectivity(x, nms=rownames(pluton))
    expect_identical(ac * totN, totConn, label="averageConnectivity: manual calculations match output")
  } # Close for (i in 1:)
})

test_that("rscale gives sensible outputs", {
  for (i in 1:100){
    x <- rt(100, df=runif(1, 1, 8)) * runif(1, 1, 100) + runif(1, 0, 10)
    xs <- rscale(x)

    expect_equal(median(xs), 0, label="rscale centers data as expected")
    expect_equal(mad(xs), 1, label="rscale scales data as expected")

    x <- x + 100
    xs <- rscale(x, center = FALSE)

    expect_gt(median(xs), 0, label="rscale doesn't center when told not to")
    expect_equal(mad(xs), 1, label="rscale uncentered output scaled as expected")

    xs <- rscale(x, scale=FALSE)

    expect_equal(median(xs), 0, label="rscale unscaled output centered as expected")
    expect_gt(mad(xs), 1, label="rscale doesn't scale when told not to")

    xs <- rscale(x, center=FALSE, scale=FALSE)
    expect_identical(x, xs, label="rscale neither centers nor scales when told not to")

    # Missing values
    nasize <- sample(3:50, size=1)
    ii <- sample(1:100, size=nasize, replace=FALSE)
    x[ii] <- NA
    xs <- rscale(x)
    expect_equal(length(xs), length(x), label="rscale returns object of same length as input when there are NAs")
    expect_equal(sum(is.na(xs)), nasize, label="rscale returns object with correct number of missing values")

    expect_error(xs <- rscale(x, na.action=na.fail), label="rscale fails with NAs when na.action=na.fail")
  }
})
