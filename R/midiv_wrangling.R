#' @name filterOTU
#' @title Filter OTUs
#'
#' @description Filtering OTUs in a midiv object.
#'
#' @param midiv.obj A midiv object (list), see \code{\link{midivObject}}.
#' @param OTU_to_keep A logical vector indicating which OTUs to keep (TRUE).
#'
#' @details This function filters OTUs, i.e. keep a subsets of rows in the
#' \code{readcount.mat} and \code{sequence.tbl} inside the \code{midiv.obj} list.
#'
#' The \code{OTU_to_keep} should be a logical vector with one element for each
#' OTU in the object.
#'
#' @return A midiv object (list) with (potentially) fewer OTUs.
#'
#' @author Lars Snipen.
#'
#' @export filterOTU
#'
filterOTU <- function(midiv.obj, OTU_to_keep){
  keep.rows <- which(OTU_to_keep)
  midiv.obj$readcount.mat <- midiv.obj$readcount.mat[keep.rows,]
  midiv.obj$sequence.tbl <- midiv.obj$sequence.tbl[keep.rows,]
  return(midiv.obj)
}


#' @name filterSample
#' @title Filter samples
#'
#' @description Filtering samples in a midiv object.
#'
#' @param midiv.obj A midiv object (list), see \code{\link{midivObject}}.
#' @param samples_to_keep A logical vector indicating which sampless to keep (TRUE).
#'
#' @details This function filters samples, i.e. keep a subsets of rows in the
#' metadata.tbl and corresponding columns in the \code{readcount.mat} of the
#' \code{midiv.obj} list.
#'
#' The \code{samples_to_keep} should be a logical vector with one element for each
#' sample in the object.
#'
#' @return A midiv object (list) with (potentially) fewer samples.
#'
#' @author Lars Snipen.
#'
#' @export filterSample
#'
filterSample <- function(midiv.obj, samples_to_keep){
  keep.cols <- which(samples_to_keep)
  midiv.obj$readcount.mat <- midiv.obj$readcount.mat[,keep.cols]
  midiv.obj$metadata.tbl <- midiv.obj$metadata.tbl[keep.cols,]
  return(midiv.obj)
}
