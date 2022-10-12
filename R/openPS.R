#' @name openPS
#' @title Creating a phyloseq-like object
#'
#' @description Creating a an object similar to a phyloseq-object, but as an open list.
#'
#' @param metadata.file Full name of text-file with metadata.
#' @param readcounts.file Full name of text-file with readcount data.
#' @param centroids.file.Full name of fasta-file with sequences.
#' @param taxonomy.file Full name of text-file with the taxonomy table results (optional).
#'
#' @details This function reads the input files, and store their content in a
#' \code{list} that we refer to as an *openPS object*.
#'
#' @return A list with the elements:
#'
#' * sample_data.tbl - a data.frame with one row for each sample. Its rownames are the SampleID's.
#' * otu_table.mat - a matrix with the readcounts. These may later be transformed. The samples are in the columns, the OTUs in the rows.
#' * sequence.tbl - a data.frame with the OTU sequences (see \code{\link{readFasta}}).
#' * tax_table.mat - a matrix with the taxonomy, one row for each OTU.
#'
#' @author Lars Snipen.
#'
#' @importFrom dplyr %>% select mutate if_else
#' @importFrom microseq readFasta
#'
#' @export openPS
#'
openPS <- function(metadata.file, readcount.file, centroids.file,
                            taxonomy.file = NULL, sintax.threshold = 0.8){
  meta.tbl <- read.table(metadata.file, header = T, sep = "\t")
  rownames(meta.tbl) <- meta.tbl$SampleID
  rc.tbl <- read.table(readcount.file, header = T, sep = "\t", comment.char = "")
  rc.mat <- rc.tbl %>%
    select(-1) %>%
    as.matrix()
  rownames(rc.mat) <- rc.tbl[,1]
  centroids.tbl <- readFasta(centroids.file)
  openPS.obj <- list(sample_data.tbl = meta.tbl,
                     otu_table.mat = rc.mat,
                     sequence.tbl = centroids.tbl)
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
    openPS.obj$tax_table.mat <- tax.mat
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
#' \code{list}, to a \code{phyloseq} object from the \code{phyloseq} R package.
#'
#' @return A \code{phyloseq} object.
#'
#' @author Lars Snipen.
#'
#' @importFrom phyloseq phyloseq otu_table sample_data tax_table
#'
#' @export openPS2phyloseq
#'
openPS2phyloseq <- function(openPS.obj){
  if(exists("tax_table.mat", where = openPS.obj)){
    if(exists("tree", where = openPS.obj)){
      ps.obj <- phyloseq(otu_table(openPS.obj$otu_table.mat, taxa_are_rows = T),
                          sample_data(openPS.obj$sample_data.tbl),
                          tax_table(openPS.obj$tax_table.mat),
                         openPS.obj$tree)
    } else {
      ps.obj <- phyloseq(otu_table(openPS.obj$otu_table.mat, taxa_are_rows = T),
                          sample_data(openPS.obj$sample_data.tbl),
                          tax_table(openPS.obj$tax_table.mat))
    }
  } else {
    if(exists("tree", where = openPS.obj)){
      ps.obj <- phyloseq(otu_table(openPS.obj$otu_table.mat, taxa_are_rows = T),
                          sample_data(openPS.obj$sample_data.tbl),
                         openPS.obj$tree)
    } else {
      ps.obj <- phyloseq(otu_table(openPS.obj$otu_table.mat, taxa_are_rows = T),
                          sample_data(openPS.obj$sample_data.tbl))
    }
  }
  return(ps.obj)
}
