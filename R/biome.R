#' @name biom2generic
#' @title Converts biom-object
#'
#' @description The content of a biom-object is converted into a readcount-table
#' and a fasta-object and written to files.
#'
#' @param biom.file Name of a biom-file output from DEBLUR,
#' @param out.folder Name of folder to output the resulting files.
#'
#' @details This function is used in our DEBLUR pipeline, to convert the
#' output from the biom-file to more generic files containing a readcount-
#' table and sequences in fasta format.
#'
#' @return The function will produce two files: The file named
#' \code{readcount_table.txt} containing the tab-separated tale of readcounts,
#' one row for each OTU and and column for each sample. The first column is
#' the names of the OTU's. The second file is named \code{centroids.fasta} and
#' contains the OTU-sequences, i.e. it has the same number of entries as there
#' are rows in the \code{readcount_table.txt}.
#'
#' @author Lars Snipen.
#'
#' @importFrom biomformat read_biom
#' @importFrom microseq writeFasta
#'
#' @export biom2generic
#'
biom2generic <- function(biom.file){
  biom.obj <- read_biom(normalizePath(biom.file))
  readcount.mat <- matrix(unlist(biom.obj$data),
                          nrow = length(biom.obj$rows),
                          ncol = length(biom.obj$columns),
                          byrow = T)
  colnames(readcount.mat) <- sapply(biom.obj$columns, function(x){x$id})
  readcount.tbl <- data.frame(DBLR = paste0("DBLR", 1:length(biom.obj$rows)))
  readcount.tbl <- cbind(readcount.tbl, as.data.frame(readcount.mat))
  write.table(readcount.tbl, file = "readcounts_deblur.txt",
              quote = F, sep = "\t", row.names = F, col.names = T)
  centroids.fsa <- data.frame(Header = paste0("DBLR", 1:length(biom.obj$rows)),
                              Sequence = sapply(biom.obj$rows, function(x){x$id}))
  writeFasta(centroids.fsa, out.file = "centroids_deblur.fasta")
}
