context("consensus")

test_that("cluster consensus behaves as expected", {
  # Test that cluster consensus is poor when there are no true clusters, and
  # good when there are true clusters. These are pretty weak tests. Doing
  # calcICL(ConsensusClusterPlus(x)) gives results with NaNs in them, which
  # doesn't seem right.
  set.seed(20170315)

  # Random uniform clustering should result in cluster consensus close to 1/k
  # for each cluster
  # Put thresh = 0.75 for k=2 and thresh = 0.25 for k=10. Implies the following
  uthresh <- 7/8 - 1/16 * 2:10

  for (i in 1:10){
    x <- matrix(runif(1000), ncol=10)
    cc <- conclus(daisy(x), K=10)

    con <- consensus(cc, type="cluster")
    expect_equal(min(c(con)), 0, label="clusterConsensus: minimum is 0 for noise")
    expect_lte(max(c(con)), 1, label="clusterConsensus: maximum is <= 1 for noise")

    con <- apply(con, 2, function(X) X / sum(X))
    # If random, each non-zero element should be close to 1/k

    con <- apply(con, 2, max)
    expect_true(all(con < uthresh), label="clusterConsensus: random uniform clustering looks random")

    # Generate some well separated cluster data
    # Create random cluster sizes
    n <- sample(10:20, size=4)
    n <- c(n, 100 - sum(n))

    # Generate 5 5-d mean vectors (the rows of mu)
    mu <- matrix(rnorm(5 * 5), ncol=5)
    mu <- mu + seq(0, 12, by=3)
    # There must be an efficient and pretty way of doing the following
    mu <- do.call("rbind", lapply(1:length(n), function(X){
      matrix(rep(mu[X, ], n[X]), ncol=ncol(mu), byrow=TRUE)
    }))
    mu <- cbind(mu, matrix(0, ncol=5, nrow=100))

    x <- matrix(rnorm(1000), ncol=10) + mu
    cc <- conclus(daisy(x), K=10)

    con <- consensus(cc, type="cluster")
    expect_equal(min(c(con)), 0, label="clusterConsensus: minimum is 0 for true clusters")
    expect_lte(max(c(con)), 1, label="clusterConsensus: maximum is <= 1 for true clusters")

    # Get variance of non-zero values. Should dive upwards after k=5
    v <- apply(con, 2, function(X) var(X[X > 0]))
    #plot(v)

    expect_gt(mean(v[5:9]), v[4] , label="clusterConsensus: variance dives upwards where there are too many clusters")
  }
})

test_that("item consensus behaves as it should", {
  data(ii)
  K <- 3 # Riskier as K gets higher because cluster IDs are arbitrary
  cc <- conclus(dist(pluton), K=K, R=1000)
  cci <- consensus(cc, type="item")

  expect_lte(max(cci), 1, label="itemConsensus: maximum value is <=1")
  expect_equal(min(cci), 0, label="itemConsensus: minimum value is 0")
  expect_equal(dim(cci), c(nrow(pluton), K, K - 1), label="itemConsensus: dimensions are as expected")

  for (k in cc$K){
    # Restructure the relevant part of ii and put it in the same order as cci
    ccp <- ii[ii[, 1] == k, -1] %>%
      mutate(itemConsensus = ifelse(is.nan(itemConsensus), 1, itemConsensus)) %>%
      mutate(item = as.numeric(as.character(item))) %>%
      spread(cluster, itemConsensus)

    for (i in 1:k){
      #print(paste(k, i, cor(cci[, i, k-1], ccp[, i+1])))
      expect_gt(cor(cci[, i, k-1], ccp[, i+1]), .9, label=paste0("itemConsensus: k=", k, ", i=", i, " - correlation high"))
    }
  }
})
