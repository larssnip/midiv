#' @name contrast_cols
#' @title Contrasting colors
#'
#' @description Making contrasting colors for use in bar charts.
#'
#' @usage contrast_cols(n_colors = 3, palette1 = "Pastel 1", palette2 = "Dark 3")
#'
#' @param n_colors The number of colors you want (integer).
#' @param palette1 A qualitative HCL palette from the \code{colorspace} package.
#' @param palette2 A qualitative HCL palette from the \code{colorspace} package.
#'
#' @details The default colors used by \code{\link{ggplot2}} are too similar
#' for neighboring sectors in bar charts, making it difficult to see where one
#' sector ends and the next starts.
#'
#' This function creates an alternative vector of colors, by using two different
#' palettes, going from lighter to darker colors, then reverse one of them and
#' finally mixing them taking every second color from each. This mean any two
#' neighboring colors are quite distinct from each other.
#'
#' You may change the colors by supplying alternative qualitative HCL palettes
#' from the \code{colorspace} package, see the vignette of this package for
#' more on what you may choose.
#'
#' @return A vector of length \code{n_colors} with colors.
#'
#' @author Lars Snipen.
#'
#' @examples
#'
#' @importFrom colorspace qualitative_hcl
#'
#' @export contrast_cols
#'
contrast_cols <- function(n_colors = 3, palette1 = "Pastel 1", palette2 = "Dark 3"){
  require(colorspace)
  p1 <- qualitative_hcl(ceiling(n_colors/2), palette = palette1)
  p2 <- qualitative_hcl(ceiling(n_colors/2), palette = palette2)
  return(as.vector(matrix(c(p1, rev(p2)), nrow = 2, byrow = T))[1:n_colors])
}
