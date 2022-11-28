#' @name nonparDA
#' @title Non-parametric DA testing
#'
#' @description Testing for differential abundance using non-parametric tests.
#'
#' @usage nonparDA(ps.obj, category_column = NULL, contrast = NULL, p.adjust.method = "BH", verbose = TRUE)
#'
#' @param ps.obj A phyloseq object.
#' @param category_column Name of one column in sample_table(ps.obj) used to group the samples.
#' @param contrast Optional specification of which category levels to use, see below.
#' @param p.adjust.method The method used for multiple testing correction, see \code{\link{p.adjust}}.
#' @param verbose Logical to turn on/off output during computing.
#'
#' @details Performs a Kruskal-Wallis non-parametric test for differential abundance
#' for each OTU in a \code{\link{phyloseq}} object.
#'
#' The \code{category_column} must be the name of a column in \code{sample_table(ps.obj)}
#' that splits the samples into groups.
#'
#' If no \code{contrast} is specified, a Kruskal-Wallis test is performed, using all category levels,
#' i.e. it tests if the abundance for at least one level deviates from at least one other level.
#' If \code{contrast} contains one text it must one of the levels in \code{category_column}, and then
#' the test is contrasting this level against all the others (A versus not A). If \code{contrast} contains
#' two texts, both must be levels in \code{category_column}, and the test is contrasting the samples from
#' these two levels only (A versus B).
#'
#' @return A table with the columns
#' \itemize{
#'   \item{OTU}
#'   \item{statistic} The Kruskal-Wallis test statistic
#'   \item{p.value} Raw p-value.
#'   \item{p.adj} Adjusted p-value due to multiple testing
#'  }
#'
#' @author Lars Snipen.
#'
#' @examples
#'
#' @export nonparDA
#'
nonparDA <- function(ps.obj, category_column = NULL, contrast = NULL,
                     p.adjust.method = "BH", verbose = TRUE){
  if(is.null(category_column)) stop("You must name a column in sample_data(ps.obj) with categories used to group samples!")
  readcount.tbl <- as.data.frame(t(otu_table(ps.obj))) %>%
    mutate(SampleID = rownames(.))
  full.tbl <- as.data.frame(as.matrix(sample_data(ps.obj))) %>%
    mutate(SampleID = rownames(.)) %>%
    select(SampleID, all_of(category_column)) %>%
    left_join(readcount.tbl, by = "SampleID")
  if(!is.null(contrast)){
    if(length(contrast) == 2){
      full.tbl <- full.tbl %>%
        filter(.data[[category_column]] %in% contrast)
      tt <- table(full.tbl[,2])
      if(length(tt) != 2) stop("Cannot find both the levels: ", contrast[1], " and ", contrast[2], " in the column ", category_column)
    } else {
      full.tbl <- full.tbl %>%
        mutate(new_cat = if_else(.data[[category_column]] == contrast, as.character(contrast), str_c("not_", contrast))) %>%
        select(-all_of(category_column)) %>%
        select(SampleID, new_cat, everything())
    }
  }
  if(verbose) cat("nonparDA:\nHave", nrow(full.tbl), "samples and", ncol(full.tbl) - 2, "OTUs\n")
  result.tbl <- tibble(OTU = colnames(full.tbl)[-c(1,2)],
                       statistic = 0,
                       p.value = -1,
                       p.adj = -1)
  for(i in 1:nrow(result.tbl)){
    tst <- kruskal.test(x = full.tbl[,2+i], g = full.tbl[,2])
    result.tbl$statistic[i] <- tst$statistic
    result.tbl$p.value[i] <- tst$p.value
  }
  result.tbl$p.adj <- p.adjust(result.tbl$p.value, method = p.adjust.method)
  return(result.tbl)
}
