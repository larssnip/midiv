#' @name min_var_down_sampling
#' @title Down-sampling read counts
#'
#' @description Down-sampling read counts with a minimum variance.
#'
#' @usage min_var_down_sampling(sample.readcounts, depth = sum(sample.readcount))
#'
#' @param sample.readcounts matrix with readcount data.
#' @param size fixed number of reads in sample after down-sampling.
#' @param rngseed Random generator seed, for exact reproducibility.
#'
#' @details Down-sampling (rarefying) of read counts means reducing the number of
#' reads in all samples to a fixed target size. This is sometimes done prior to
#' some downstream analyses of read count data.
#'
#' Such down-sampling is typically done by sampling either with or without
#' replacement. The most common is with replacement, being faster. This means the
#' fixed number of reads are sampled with the original read count relative
#' frequencies as probabilities for the various outcomes. In such cases the
#' OTU's with very low read counts originally may come with more reads after
#' down-sampling.
#' The sampling without replacement means the actual reads are sampled, and no
#' OTU can end up with more reads after down-sampling. However, both procedures
#' introduce an extra variance in the data, especially in the more abundant OTU's.
#' This variance has no biological information.
#'
#' This function minimizes implements a down-sampling that minimizes the
#' variance, by only sampling the few reads needed to 'correct' the expected read
#' count into an actual read count, as follows:
#'
#' We first compute the expected read count given target \code{size}. This is
#' \code{E = depth * sample.readcounts/sum(sample.readcounts)}. From this we get
#' the base count \code{b = floor(E)}. The remainder is
#' \code{r = E - b}. The remaining reads
#' \code{depth - sum(base)} are finally distributed over
#' the OTU's using \code{r/sum(r)} as the probabilities. Thus, only the
#' (few) remaining reads will vary randomly between repeated use of this function
#' on the same data.
#'
#' @return A vector of same size as the input, but with down-sampled readcounts.
#'
#' @author Lars Snipen.
#'
#' @examples
#'
#' @export min_var_down_sampling
#'
min_var_down_sampling <- function(sample.readcounts, size = sum(sample.readcounts),
                                  rngseed = NULL){
  E <- size * sample.readcounts / sum(sample.readcounts)
  b <- floor(E)
  r <- E - b
  if(!is.null(rngseed)){
    set.seed(rngseed)
  }
  if(sum(r) > 0){
    cfac <- factor(1:length(sample.readcounts))
    rest <- table(sample(cfac, size = size - sum(b), replace = T, prob = r))
  } else {
    rest <- r
  }
  sample.downsampled <- b + rest
  names(sample.downsampled) <- names(sample.readcounts)
  return(sample.downsampled)
}

#' @name rarefaction_curves
#' @title Computes rarefaction curves
#'
#' @description Computes rarefaction curves for a phyloseq object.
#'
#' @usage rarefaction_curves(phyloseq.object, step = 100, plot = TRUE)
#'
#' @param phyloseq.object A phyloseq object.
#' @param step The step-size used to choose different down-sampling levels (readcounts).
#' @param plot Logical, plotting the curves or return the table for later plotting.
#'
#' @details This function takes as input a phyloseq object, and uses the function
#' \code{\link{rarecurve}} in the \code{vegan} package for the computations.
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
#' @importFrom dplyr rename
#'
#' @export rarefaction_curves
#'
rarefaction_curves <- function(phyloseq.object, step = 100, plot = TRUE){
  otu.tbl <- as.data.frame(otu_table(phyloseq.object))
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



