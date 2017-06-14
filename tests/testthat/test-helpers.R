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

test_that("trim.hilo behaves as expected", {
  x <- rnorm(100)
  expect_equal(length(trim.hilo(x, trim=.05)), length(x) * (1 - 2 * .05), label="trim.hilo removes 2*5% correctly")
  expect_equal(length(trim.hilo(x, trim=.15)), length(x) * (1 - 2 * .15), label="trim.hilo removes 2*15% correctly")

  for (i in 1:100){
    x <- rnorm(round(runif(1, 30, 3000)))
    a <- runif(1, 0, .4)

    res <- (length(x) * (1 - 2*a)) - length(trim.hilo(x, a))
    expect_lt(res, 2, label="trim.hilo removes random amounts correctly")
  }
})

test_that("rskew gives sensible outputs", {
  for (i in 1:10){ # this tests 95% confidence intervals so will fail occasionally
    fun <- function(){
      x <- rnorm(1000, 10)
      # Next line should give 0, > 0, < 0
      if (any(is.na(log(x)))) browser()
      c(rskew(x, trim=0), rskew(exp(x), trim=0), rskew(log(x), trim=0))
    }

    rr <- t(replicate(1000, fun()))
    ci <- apply(rr, 2, quantile, prob=c(.025, .975))

    expect_true(ci[1, 1] < 0 & ci[2, 1] > 0, label="rskew: Symmetric data have zero skew")
    expect_true(min(ci[, 2]) > 0, label="rskew: Positively skewed data have positive skew")
    expect_true(max(ci[, 3]) < 0, label="rskew: Negatively skewed data have negative skew")
  }
})

test_that("rscale gives sensible outputs", {
  for (i in 1:100){
    x <- rt(100, df=runif(1, 1, 8)) * runif(1, 1, 100) + runif(1, 0, 10)
    x <- rscale(x)

    expect_lt(median(x), 10^{-6}, label="rscale centers data as expected")
    expect_lt(diff(range(trim.hilo(x, trim=.1))), 6, label="rscale scales data as expected")
  }
})

test_that("deskew gives sensible outputs", {
  for (i in 1:100){
    x <- rnorm(1000, 4)
    a <- sample(c(-2, -1, -.5, 0.001, .5, 2), size=1)
    b <- runif(1, 1, 10)

    y <- ((x + b)^a - 1) / a

    dy <- deskew(y)
    #abs(rskew(dy, trim=0)) - abs(rskew(y, trim=0))
    expect_lt(abs(rskew(dy, trim=0)) - abs(rskew(y, trim=0)), .2, label="deskew reduces skew")
  }
})

