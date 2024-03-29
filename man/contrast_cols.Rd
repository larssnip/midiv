% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colors.R
\name{contrast_cols}
\alias{contrast_cols}
\title{Contrasting colors}
\usage{
contrast_cols(n_colors = 3, palette1 = "Pastel 1", palette2 = "Dark 3")
}
\arguments{
\item{n_colors}{The number of colors you want (integer).}

\item{palette1}{A qualitative HCL palette from the \code{colorspace} package.}

\item{palette2}{A qualitative HCL palette from the \code{colorspace} package.}
}
\value{
A vector of length \code{n_colors} with colors.
}
\description{
Making contrasting colors for use in bar charts.
}
\details{
The default colors used by \code{\link{ggplot2}} are too similar
for neighboring sectors in bar charts, making it difficult to see where one
sector ends and the next starts.

This function creates an alternative vector of colors, by using two different
palettes, going from lighter to darker colors, then reverse one of them and
finally mixing them taking every second color from each. This mean any two
neighboring colors are quite distinct from each other.

You may change the colors by supplying alternative qualitative HCL palettes
from the \code{colorspace} package, see the vignette of this package for
more on what you may choose.
}
\author{
Lars Snipen.
}
