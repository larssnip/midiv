#' @name rarefaction_curves
#' @title Computes rarefaction curves
#'
#' @description Computes rarefaction curves for a phyloseq object.
#'
#' @usage rarefaction_curves(phyloseq.object, step = 100, plot = TRUE)
#'
#' @param midiv.object A midiv object. You may supply either this or a phyloseq object.
#' @param phyloseq.object A phyloseq object. You may supply either this or a midiv object.
#' @param step The step-size used to choose different down-sampling levels (readcounts).
#' @param plot Logical, plotting the curves or return the table for later plotting.
#'
#' @details This function takes as input either a midiv object or a phyloseq object,
#' and uses the function \code{\link{rarecurve}} in the \code{vegan} package for
#' the computations.
#'
#' @return If \code{plot = TRUE} a \code{ggplot} object, if \code{plot = FALSE}
#' the table with the data for plotting the curves. The latter is useful if you want
#' a different display of the curves
#'
#' @author Lars Snipen.
#'
#' @examples
#'
#' @importFrom vegan rarecurve
#' @importFrom ggplot2 ggplot geom_line labs
#' @importFrom dplyr rename %>%
#'
#' @export rarefaction_curves
#'
rarefaction_curves <- function(midiv.object, phyloseq.object = NULL, step = 100, plot = TRUE){
  if(is.null(phyloseq.object)){
    otu.tbl <- midiv.object$readcount.mat
  } else {
    otu.tbl <- as.data.frame(otu_table(phyloseq.object))
  }
  rar.tbl <- rarecurve(t(otu.tbl), step = step, tidy = T) %>%
    rename(Sample = Site, Reads = Sample, OTUs = Species)
  if(plot){
    fig <- ggplot(rar.tbl) +
      geom_line(aes(x = Reads, y = OTUs, color = Sample)) +
      labs(x = "Number of reads", y = "Number of OTUs")
    print(fig)
    return(fig)
  } else {
    return(rar.tbl)
  }
}



