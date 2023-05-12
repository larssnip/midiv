#' @name midivObject
#' @title Creating a midiv-object
#'
#' @description Creating a list containing data from microbial community profiling at MiDiv.
#'
#' @param metadata.tbl A data.frame with the sample metadata, or the name of a file.
#' @param readcount.mat A matrix with readcounts data, or the name of a file.
#' @param sequence.tbl A data.frame with the OTU sequences, or the name of a file.
#' @param taxonomy.tbl A data.frame with the OTU taxonomy table results, or the name of a file (optional).
#' @param sample_id_column Text with the name of the metadata.tbl column name that identifies samples.
#' @param filter_samples Text indicating if the metadata or the readcount table should decide the samples to keep.
#' @param filter_OTUs Text indicating if the readcount or sequence table should decide the OTUs to keep.
#'
#' @details This function stores the data structure from the processing of microbial community
#' sequencing data in a \code{list} that we refer to as a *midiv-object*.
#'
#' The first four arguments are either data structure already read into R, or names of the files to read.
#'
#' The \code{metadata.tbl} is a data.frame with one row for each sample, containing metadata
#' for each sample in the columns. It must have one column that uniquely identifies each sample.
#' This is specified in \code{sample_id_column} and is by default \code{"SampleID"}.
#' If a file name is specified, it is assumed to be a tab-delimited text file.
#'
#' The \code{readcount.mat} is a *matrix* with readcounts, one row for each OTU and
#' one column for each sample. If a file name is specified, it is assumed to be
#' a tab-delimited text file, with OTUs in the rows and the samples in the columns,
#' but where the first column contains the OTU identifying texts (typically OTU1, OTU2,...).
#'
#' The \code{sequence.tbl} is a Fasta-table with centroid sequences, see \code{\link{readFasta}},
#' or the name of a FASTA-file. These sequences are the centroid sequences for
#' each OTU, and the texts in the Header column must match the texts identifying the
#' OTUs in the \code{readcount.mat} above (typically OTU1, OTU2,...).
#'
#' The \code{taxonomy.tbl} may be supplied, and must then be a table where the first
#' column is named \code{OTU} and contains the texts identifying the OTUs (typically OTU1, OTU2,...).
#' The remaining columns should list the taxonomy at various ranks, and nothing more. Columns of scores
#' etc must been selected out of this table, see \code{\link{sintaxFilter}}.
#' If \code{taxonomy.tbl} is a file name
#' the file must be a tab-delimited text file with the columns as described above.
#' NB! In the created \code{midiv} object this is merged with the \code{sequence.tbl}.
#'
#' The argument \code{filter_samples} is only used if the samples in the \code{metadata.tbl}
#' and \code{readcount.mat} are not the same. If \code{filter_sample = "metadata"} the samples
#' in this table are kept, and the \code{readcount.mat} is trimmed accordingly.
#'
#' The argument \code{filter_OTUs} is only used if the OTUs in the \code{readcount.mat} and
#' \code{sequence.tbl} are not the same. If \code{filter_OTUs = "sequence"} the samples
#' in this table are kept, and the \code{readcount.mat} is trimmed/extended accordingly.
#'
#' @return A list with the elements:
#'
#' \itemize{
#'   \item{\code{metadata.tbl}}{ a data.frame with one row for each sample.}
#'   \item{\code{readcount.mat}}{ a matrix with the readcounts, the samples are in the columns, the OTUs in the rows.}
#'   \item{\code{sequence.tbl}}{ a data.frame with the OTU sequences (see \code{\link{readFasta}}) and taxonomy, if supplied.}
#' }
#'
#' @author Lars Snipen.
#'
#' @importFrom readr read_delim
#' @importFrom microseq readFasta
#'
#' @export midivObject
#'
midivObject <- function(metadata.tbl, readcount.mat, sequence.tbl,
                        taxonomy.tbl = NULL, sample_id_column = "SampleID",
                        filter_samples = "metadata", filter_OTUs = "sequence"){
  ## The sample metadata
  if(is.character(metadata.tbl)){
    metadata.tbl <- suppressMessages(read_delim(metadata.tbl))
  }
  metadata.tbl <- as.data.frame(metadata.tbl)
  if(!(sample_id_column %in% colnames(metadata.tbl))) stop("Cannot find column ", sample_id_column, " in metadata.tbl")

  ## The readcount matrix
  if(is.character(readcount.mat)){
    readcount.mat <- as.data.frame(suppressMessages(read_delim(readcount.mat)))
    otu.names <- readcount.mat[,1]
    readcount.mat <- as.matrix(readcount.mat[,-1])
    rownames(readcount.mat) <- otu.names
  }
  readcount.mat <- as.matrix(readcount.mat)

  ## Sample filtering
  if(filter_samples == "metadata"){
    sample_id <- metadata.tbl[,sample_id_column]
    keep.cols <- which(colnames(readcount.mat) %in% sample_id)
    readcount.mat <- readcount.mat[,keep.cols]
  } else {
    sample_id <- colnames(readcount.mat)
    keep.rows <- which(metadata.tbl[,sample_id_column] %in% sample_id)
    metadata.tbl <- metadata.tbl[keep.rows,]
  }

  ## The sequences
  if(is.character(sequence.tbl)){
    sequence.tbl <- readFasta(sequence.tbl)
  }
  colnames(sequence.tbl) <- c("OTU", "Sequence")

  ## The taxonomy
  if(!is.null(taxonomy.tbl)){
    if(is.character(taxonomy.tbl)){
      taxonomy.tbl <- suppressMessages(read_delim(taxonomy.tbl))
    }
    sequence.tbl <- full_join(sequence.tbl, taxonomy.tbl, by = "OTU")
  }

  ## Filter OTUs
  if(filter_OTUs == "sequence"){
    sequence.tbl <- drop_na(sequence.tbl, Sequence)
    keep.rows <- which(rownames(readcount.mat) %in% sequence.tbl$OTU)
    readcount.mat <- readcount.mat[keep.rows,]
    extend.rows <- which(!(sequence.tbl$OTU %in% rownames(readcount.mat)))
    if(length(extend.rows) > 0){
      extend.mat <- matrix(0, nrow = length(extend.rows), ncol = ncol(readcount.mat))
      rownames(extend.mat) <- sequence.tbl$OTU[extend.rows]
      readcount.mat <- rbind(readcount.mat, extend.mat)
    }
  }

  ## Ensure same ordering of samples
  idx <- match(metadata.tbl[,sample_id_column], colnames(readcount.mat))
  if(sum(is.na(idx)) > 0) stop("The SampleID's in ", metadata.file, " do not correspond to those in ", readcounts.file)
  readcount.mat <- readcount.mat[,idx]
  ## Ensure same ordering of OTUs
  idx <- match(sequence.tbl$OTU, rownames(readcount.mat))
  if(sum(is.na(idx)) > 0) stop("The OTU's in ", centroids.file, " do not correspond to those in ", readcounts.file)
  readcount.mat <- readcount.mat[idx,]

  ## Create list
  midiv.obj <- list(metadata.tbl = metadata.tbl,
                     readcount.mat = readcount.mat,
                     sequence.tbl = sequence.tbl)
  return(midiv.obj)
}



#' @name midiv2phyloseq
#' @title Convert midiv to phyloseq object
#'
#' @description Creating a phyloseq object from a midiv object.
#'
#' @param midiv.obj A midiv object, see \code{\link{midivObject}}.
#' @param sample_id_column Text with the name of the metadata.tbl column name that identifies samples.
#'
#' @details This function converts a midiv object, which is a simple
#' \code{list}, to a \code{\link{phyloseq}} object from the \code{phyloseq} R package.
#'
#' @return A \code{\link{phyloseq}} object.
#'
#' @author Lars Snipen.
#'
#' @importFrom phyloseq phyloseq otu_table sample_data tax_table
#' @importFrom dplyr %>% select mutate
#'
#' @export midiv2phyloseq
#'
midiv2phyloseq <- function(midiv.obj, sample_id_column = "SampleID"){
  otu.table <- midiv.obj$readcount.mat

  sample.data <- midiv.obj$metadata.tbl
  rownames(sample.data) <- midiv.obj$metadata.tbl[,sample_id_column]

  taxonomy.tbl <- select(midiv.obj$sequence.tbl, -c(OTU, Sequence))
  if(ncol(taxonomy.tbl) > 0){
    tax.mat <- as.matrix(taxonomy.tbl)
    rownames(tax.mat) <- midiv.obj$sequence.tbl$OTU
    ps.obj <- phyloseq(otu_table(otu.table, taxa_are_rows = T),
                       sample_data(sample.data),
                       tax_table(tax.mat))
  } else {
    ps.obj <- phyloseq(otu_table(otu.table, taxa_are_rows = T),
                       sample_data(sample.data))
  }
  return(ps.obj)
}


#' @name phyloseq2midiv
#' @title Convert phyloseq to midiv object
#'
#' @description Creating a simple list from a phyloseq object.
#'
#' @param phyloseq.obj A phyloseq object, see \code{\link{phyloseq}}.
#'
#' @details This function converts a phyloseq object to a simple \code{\link{list}}
#' giving the entries the names as in \code{\link{midivObject}}.
#'
#' This may be convenient for some wrangling on the data, and then perhaps converting
#' it back to a phyloseq object again with \code{\link{midiv2phyloseq}}.
#'
#' @return A \code{list} with entries as in a midiv object, except that the
#' \code{sequence.tbl} do not contain sequences, only taxonomy.
#'
#' @author Lars Snipen.
#'
#' @importFrom phyloseq phyloseq otu_table sample_data tax_table
#'
#' @export phyloseq2midiv
#'
phyloseq2midiv <- function(phyloseq.obj){
  lst <- list(metadata.tbl = as.data.frame(as.matrix(sample_data(phyloseq.obj))),
              readcount.mat = as.matrix(as.data.frame(otu_table(phyloseq.obj))),
              sequence.tbl = as.data.frame(tax_table(phyloseq.obj)))
  return(lst)
}
