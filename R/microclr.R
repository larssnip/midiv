#' @name clr
#' @title Centred log-ratio transform
#'
#' @description Transforms readcount data with Aitchisons transform.
#'
#' @usage clr(sample.readcounts, n.pseudo = 1)
#'
#' @param sample.readcount A vector with readcount data for one sample.
#' @param n.pseudo number of pseudo-readcounts to add.
#'
#' @details This is a standard implementation of the Aitchisons centered log-ratio
#' transform (Aitchison 1986) for compositional data. Readcount data can be seen as
#' compositional data since the total number of readcounts in a sample does not
#' carry any information about the biology, but is simply an effect of sequencing
#' depth. Thus, the information in the data lies in the relative values, not the
#' absolute. By transforming such data with this function, you get data who are
#' better suited for a number of downstream analyses, e.g. typically analyses
#' making use of sum-of-squares type of statistics, like PCA, PLS, ANOVA or
#' clustering with euclidean distances.
#'
#' The transform does not accept zeros in any cell of the \code{sample.readcounts}. To
#' cope with this you add pseudo-counts. By default 1 additional readcount is assigned to
#' all cells in \code{readcount.mat}. You may change this value, and it need not
#' be an integer. The rationale behind this is that we *a priori* assumes a uniform
#' distribution of the taxa, and the more pseudo-counts you add, the more weight
#' you give to this prior.
#'
#' The main reason this takes a vector as input, not a full matrix, is that
#' readcount data are stored in several ways, bot in tables and matrices. Also,
#' having a per-sample oriented implementation makes it straightforward to use
#' this in conjunction with the \code{phyloseq} R package, and the function
#' \code{\link{transform_sample_counts}}
#'
#' @return A vector of same size as the input, but with transformed readcounts.
#'
#' @author Lars Snipen.
#'
#' @references Aitchison J. The Statistical Analysis of Compositional Data. London, UK: Chapman & Hall; 1986.
#'
#' @examples
#'
#' @export clr
#'
clr <- function(sample.readcounts, n.pseudo = 1){
  sample.readcounts <- sample.readcounts + n.pseudo
  sample.rel.abd <- sample.readcounts / sum(sample.readcounts)
  lc <- log2(sample.rel.abd)
  clr.values <- lc - mean(lc)
  names(clr.values) <- names(sample.readcounts)
  return(clr.values)
}

