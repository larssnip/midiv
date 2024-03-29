% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/midiv_object.R
\name{midiv2phyloseq}
\alias{midiv2phyloseq}
\title{Convert midiv to phyloseq object}
\usage{
midiv2phyloseq(midiv.obj, sample_id_column = "SampleID")
}
\arguments{
\item{midiv.obj}{A midiv object, see \code{\link{midivObject}}.}

\item{sample_id_column}{Text with the name of the metadata.tbl column name that identifies samples.}
}
\value{
A \code{\link{phyloseq}} object.
}
\description{
Creating a phyloseq object from a midiv object.
}
\details{
This function converts a midiv object, which is a simple
\code{list}, to a \code{\link{phyloseq}} object from the \code{phyloseq} R package.
}
\author{
Lars Snipen.
}
