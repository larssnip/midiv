#' @name openPS
#' @title Creating a phyloseq-like object
#'
#' @description Creating an object similar to a phyloseq-object, but as an open list.
#'
#' @param metadata.file Full name of text-file with sample metadata. This is named sample_table in phyloseq.
#' @param readcounts.file Full name of text-file with readcounts data, assume samples in the columns. This is named otu_table (but is a matrix) in phyloseq.
#' @param centroids.file Full name of fasta-file with OTU sequences. This is not found in phyloseq.
#' @param taxonomy.file Full name of text-file with the OTU taxonomy table results (optional). This is named tax_table (but is a matrix) in phyloseq.
#' @param sintax.threshold A value between 0.0 and 1.0 indicating confidence threshold for the taxonomy (optional).
#'
#' @details This function reads the input files, and store their content in a
#' \code{list} that we refer to as an *openPS object*.
#'
#' If a \code{taxonomy.file} is specified, the \code{sintax.threshold} will be used to assign
#' all classification with a score below this threshold to \code{"unclassified"} in the \code{tax_table.mat}.
#'
#' @return A list with the elements:
#'
#' \itemize{
#'   \item{\code{metadata.tbl}}{ a data.frame with one row for each sample. Its rownames are the SampleID's.}
#'   \item{\code{readcount.mat}}{ a matrix with the readcounts. These may later be transformed. The samples are in the columns, the OTUs in the rows.}
#'   \item{\code{sequence.tbl}}{ a data.frame with the OTU sequences (see \code{\link{readFasta}}).}
#'   \item{\code{taxonomy.mat}}{ a matrix with the taxonomy, one row for each OTU.}
#' }
#'
#' @author Lars Snipen.
#'
#' @importFrom readr read_delim
#' @importFrom dplyr %>% select mutate if_else
#' @importFrom microseq readFasta
#'
#' @export openPS
#'
openPS <- function(metadata.file, readcounts.file, centroids.file,
                   taxonomy.file = NULL, sintax.threshold = 0.0){
  ## The sample metadata
  metadata.tbl <- suppressMessages(read_delim(metadata.file, delim = "\t")) %>%
    as.data.frame()
  rownames(metadata.tbl) <- metadata.tbl$SampleID  # convenient for converting to phyloseq

  ## The readcounts as a matrix
  readcount.tbl <- suppressMessages(read_delim(readcounts.file, delim = "\t")) %>%
    as.data.frame()
  readcount.mat <- readcount.tbl %>%
    select(-1) %>%
    as.matrix()
  rownames(readcount.mat) <- readcount.tbl[,1]

  ## The sequences
  sequence.tbl <- readFasta(centroids.file)

  ## Ensure same ordering
  idx <- match(rownames(metadata.tbl), colnames(readcount.mat))
  if(sum(is.na(idx)) > 0) stop("The SampleID's in ", metadata.file, " do not correspond to those in ", readcounts.file)
  readcount.mat <- readcount.mat[,idx]
  idx <- match(sequence.tbl$Header, rownames(readcount.mat))
  if(sum(is.na(idx)) > 0) stop("The OTU's in ", centroids.file, " do not correspond to those in ", readcounts.file)
  readcount.mat <- readcount.mat[idx,]

  ## Create list
  openPS.obj <- list(metadata.tbl = metadata.tbl,
                     readcount.mat = readcount.mat,
                     sequence.tbl = sequence.tbl)
  if(!is.null(taxonomy.file)){
    tax.tbl <- read.table(taxonomy.file, header = T, sep = "\t") %>%
      mutate(species = as.character(species)) %>%
      mutate(species = if_else(species_score >= sintax.threshold, species, "unclassified")) %>%
      mutate(genus = if_else(genus_score >= sintax.threshold, genus, "unclassified")) %>%
      mutate(family = if_else(family_score >= sintax.threshold, family, "unclassified")) %>%
      mutate(order = if_else(order_score >= sintax.threshold, order, "unclassified")) %>%
      mutate(class = if_else(class_score >= sintax.threshold, class, "unclassified")) %>%
      mutate(phylum = if_else(phylum_score >= sintax.threshold, phylum, "unclassified")) %>%
      mutate(domain = if_else(domain_score >= sintax.threshold, domain, "unclassified")) %>%
      select(-contains("score"))
    tax.mat <- tax.tbl %>%
      select(-1) %>%
      as.matrix()
    rownames(tax.mat) <- tax.tbl[,1]

    ## Ensure same ordering
    idx <- match(rownames(readcount.mat), rownames(tax.mat))
    if(sum(is.na(idx)) > 0) stop("The OTU's in ", taxonomy.file, "do not correspond to those in ", readcounts.file)
    tax.mat <- tax.mat[idx,]

    openPS.obj$taxonomy.mat <- tax.mat
  }
  return(openPS.obj)
}

#' @name openPS2phyloseq
#' @title Convert to phyloseq object
#'
#' @description Creating a phyloseq object from an openPS object.
#'
#' @param openPS.obj An openPS object, see \code{\link{openPS}}.
#'
#' @details This function converts an openPS object, which is a simple
#' \code{list}, to a \code{\link{phyloseq}} object from the \code{phyloseq} R package.
#'
#' @return A \code{\link{phyloseq}} object.
#'
#' @author Lars Snipen.
#'
#' @importFrom phyloseq phyloseq otu_table sample_data tax_table
#'
#' @export openPS2phyloseq
#'
openPS2phyloseq <- function(openPS.obj){
  if(exists("taxonomy.mat", where = openPS.obj)){
    if(exists("tree", where = openPS.obj)){
      ps.obj <- phyloseq(otu_table(openPS.obj$readcount.mat, taxa_are_rows = T),
                         sample_data(openPS.obj$metadata.tbl),
                         tax_table(openPS.obj$taxonomy.mat),
                         openPS.obj$tree)
    } else {
      ps.obj <- phyloseq(otu_table(openPS.obj$readcount.mat, taxa_are_rows = T),
                         sample_data(openPS.obj$metadata.tbl),
                         tax_table(openPS.obj$taxonomy.mat))
    }
  } else {
    if(exists("tree", where = openPS.obj)){
      ps.obj <- phyloseq(otu_table(openPS.obj$readcount.mat, taxa_are_rows = T),
                         sample_data(openPS.obj$metadata.tbl),
                         openPS.obj$tree)
    } else {
      ps.obj <- phyloseq(otu_table(openPS.obj$readcount.mat, taxa_are_rows = T),
                         sample_data(openPS.obj$metadata.tbl))
    }
  }
  return(ps.obj)
}


#' @name phyloseq2openPS
#' @title Convert to openPS object
#'
#' @description Creating a simple list from a phyloseq object.
#'
#' @param phyloseq.obj A phyloseq object, see \code{\link{phyloseq}}.
#'
#' @details This function converts a phyloseq object to a simple \code{\link{list}}
#' giving the entries the names as in \code{\link{openPS}}.
#'
#' This may be convenient for some wrangling on the data, and then perhaps converting
#' it back to a phyloseq object again with \code{\link{openPS2phyloseq}}.
#'
#' @return A \code{list} with entries as in an openPS object.
#'
#' @author Lars Snipen.
#'
#' @importFrom phyloseq phyloseq otu_table sample_data tax_table
#'
#' @export phyloseq2openPS
#'
phyloseq2openPS <- function(phyloseq.obj){
  lst <- list(metadata.tbl = as.data.frame(as.matrix(sample_data(phyloseq.obj))),
              readcount.mat = as.matrix(as.data.frame(otu_table(phyloseq.obj))),
              taxonomy.mat = as.matrix(as.data.frame(tax_table(phyloseq.obj))))
  return(lst)
}
