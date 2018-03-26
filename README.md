# Consensus Clustering in R #

This is the consensus clustering approach of Monti, Tamayo,
Mesirov and Golub. Features include the following.

* Users can provide their own clustering functions into conclus.
* The returned object has print, summary and ggplot methods.
* The object returned by summary has its own print and ggplot methods.
* Designed to work specifically with dissimilarity or dist matrices.
* Test suite.

**AstraZeneca** is partially funding this work.

## Updates
* 2018-03-26. If a cluster had ties in terms of representative items, the representatives function returned a list. The function now arbitrarily returns the first representative item per cluster.

## References ##
S. Monti, P. Tamayo, J. Mesirove and T. Golub, _Consensus
clustering: a resampling-based method for class discovery
and visualization of gene expression microarray data_,
Machine Learning, 52, 91 -- 118, 2003

