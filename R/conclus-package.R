#' @importFrom cluster pam fanny daisy
#' @importFrom stats hclust as.dist runif
#' @importFrom ggplot2 ggplot aes geom_tile theme facet_wrap geom_line element_blank
#'   scale_x_continuous scale_y_continuous geom_point scale_fill_gradient
#'   scale_y_discrete scale_x_discrete
#' @importFrom tidyr gather
#' @importFrom dplyr mutate bind_rows
#' @importFrom gridExtra grid.arrange
#' @importFrom graphics lines plot

globalVariables(c("%>%", "M", "vars", "vars2", "k", "CDF", "index", "pluton"))

#' St Jude leukemia data
#' @name stjude
#' @docType data
#' @format A data frame containing data 983 features (columns) measured on 246 patients (rows).
#' @source http://portals.broadinstitute.org/cgi-bin/cancer/publications/view/87 Beware! the URL
#'   has been known to change. An alternative is to go to http://portals.broadinstitute.org/cgi-bin/cancer
#'   and search for the paper by Monti et al.
#' @note It is reported that there are 6 distinct groups within the data, so something
#'   close to 6 ought to be identified by clustering either the rows or the columns.
#' @references S. Monti, P. Tamayo, J. Mesirov and T. Golub, Consensus clustering: a
#'   resampling-based method for class discovery and visualization of gene expression data,
#'   Machine Learning, 52, 91 -- 118, 2003
NULL

#' Gaussian3 simulated data
#' @name Gaussian3
#' @aliases Uniform1
#' @docType data
#' @format Data frames of simulated data, described in Monti et al.
#' @source http://portals.broadinstitute.org/cgi-bin/cancer/publications/view/87 Beware! the URL
#'   has been known to change. An alternative is to go to http://portals.broadinstitute.org/cgi-bin/cancer
#'   and search for the paper by Monti et al.
#' @references S. Monti, P. Tamayo, J. Mesirov and T. Golub, Consensus clustering: a
#'   resampling-based method for class discovery and visualization of gene expression data,
#'   Machine Learning, 52, 91 -- 118, 2003
NULL

# The following code was used to create item consensus info using ConcensusClusterPlus
if (FALSE){
  library(ConsensusClusterPlus)
  ccp <- ConsensusClusterPlus(dist(pluton), maxK=10, pItem=.5, clusterAlg="pam", reps=1000)
  ii <- calcICL(ccp)[[2]]
  save(ii, file="data/ii.rda")
}

#' Item consensus estimates from an alternative implementation
#' @name ii
#' @docType data
#' @aliases ii
#' @format A data frame containing the output of running
#'  \code{calcICL(ConsensusClusterPlus(dist(pluton), maxK=10, pItem=.5, clusterAlg="pam", reps=1000))[[2]]}
#'  only included for testing \code{itemConsensus}.
#' @source Output from \code{ConcensusClusterPlus::calcICL}
#' @keywords datasets
#' @references Wilkerson, M.D., Hayes, D.N. (2010). ConsensusClusterPlus: a class discovery tool with confidence assessments and item tracking. Bioinformatics, 2010 Jun 15;26(12):1572-3
NULL
