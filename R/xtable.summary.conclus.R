#' Create an xtable from a summary.conclus object
#' @param x The output of \code{summary.conclus}.
#' @param caption, label, align, digits, display, auto, ... Arguments to \code{xtable}
#' @method xtable summary.conclus
#' @export xtable.summary.conclus
xtable.summary.conclus <- function(x, caption=NULL, label=NULL, align=NULL, digits=2,
                                   display=NULL, auto=FALSE, ...){
  ccc <- round(x$Consensus, digits=digits)

  ccc <- matrix(as.character(ccc), ncol=ncol(ccc))
  ccc[ccc == "0"] <- ""


  dimnames(ccc) <- dimnames(x$Consensus)

  xtable(ccc, caption=caption, label=label, align=align, digits=digits, display=display, auto=auto, ...)
}
