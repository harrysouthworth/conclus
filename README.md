# Consensus Clustering in R #

This is the consensus clustering approach of Monti, Tamayo,
Mesirov and Golub. Features (will) include the following.

* Users can provide their own clustering functions into conclus.
* The returned object has print, summary and ggplot methods.
* The object returned by summary has its own print and ggplot methods.
* Preprocessing functions for (robustly) deskewing and scaling data are provided.
* Designed to work specifically with dissimilarity or dist matrices.
* Test suite under construction.

**AstraZeneca** is partially funding this work.


## References ##
S. Monti, P. Tamayo, J. Mesirove and T. Golub, _Consensus
clustering: a resampling-based method for class discovery
and visualization of gene expression microarray data_,
Machine Learning, 52, 91 -- 118, 2003

