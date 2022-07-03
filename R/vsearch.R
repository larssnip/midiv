#' @name vsearch_update_metadata
#' @title Updating metadata file
#'
#' @description Updating the metadata file with information on VSEARCH processing.
#'
#' @param metadata.file Full name of file with metadata (text file).
#' @param tmp.folder	Name of folder where fastq files are found.
#' @param readcounts.file.OTU Name of file with OTU readcounts table.
#' @param readcounts.file.ZOTU Name of file with ZOTU readcounts table.
#'
#' @details The \code{metadata.file} must be a text file with a table
#' with one row for each sample. It must follow the MiDiv metadata table standard
#' format. The columns used by this function are:
#' * ProjectID
#' * SequencingRunID
#' * SampleID
#'
#' The names of the files in \code{tmp.folder} should all follow
#' the format: ProjectID_SequencingRunID_SampleID.fq.
#'
#' This function will read the files and update the metadata file with the columns
#' n_vsearch_merged, n_vsearch_reads_OTU and n_vsearch_reads_ZOTU.
#'
#' @return The function will update the content of the file \code{metadata.file}.
#'
#' @author Lars Snipen.
#'
#' @importFrom readr read_delim write_delim
#' @importFrom dplyr %>% select mutate left_join
#' @importFrom microseq readFastq writeFastq
#' @importFrom magrittr set_colnames
#'
#' @export vsearch_update_metadata
#'
vsearch_update_metadata <- function(metadata.file, tmp.folder = "tmp_vsearch",
                                    readcounts.file.OTU = "readcounts_vsearch_OTU.txt",
                                    readcounts.file.ZOTU = "readcounts_vsearch_ZOTU.txt"){
  cat("Reading the metadata file ", metadata.file, "...\n")
  metadata.tbl <- suppressMessages(read_delim(metadata.file, delim = "\t")) %>%
    mutate(n_vsearch_merged = 0)
  if(exists("n_vsearch_reads_OTU", metadata.tbl)) metadata.tbl <- select(metadata.tbl, -n_vsearch_reads_OTU)
  if(exists("n_vsearch_reads_ZOTU", metadata.tbl)) metadata.tbl <- select(metadata.tbl, -n_vsearch_reads_ZOTU)
  cat("Looping over samples to collect number of merged reads")
  for(i in 1:nrow(metadata.tbl)){
    file.stem <- str_c(metadata.tbl$ProjectID[i], "_",
                       metadata.tbl$SequencingRunID[i], "_",
                       metadata.tbl$SampleID[i])
    fq <- readFastq(file.path(tmp.folder, str_c(file.stem, ".fq")))
    metadata.tbl$n_vsearch_merged[i] <- nrow(fq)
    cat(".")
  }
  cat("\nReading readcounts from ", readcounts.file.OTU, "...\n")
  rc.tbl <- suppressMessages(read_delim(readcounts.file.OTU, delim = "\t")) %>%
    mutate(SampleID = as.character(SampleID))
  readcount.mat <- rc.tbl %>%
    select(-1) %>%
    as.matrix() %>%
    t() %>%
    set_colnames(rc.tbl$`#OTU ID`)
  metadata.tbl <- left_join(metadata.tbl,
                            data.frame(n_vsearch_reads_OTU = rowSums(readcount.mat),
                                       SampleID = rownames(readcount.mat)),
                            by = "SampleID")
  cat("Reading readcounts from ", readcounts.file.ZOTU, "...\n")
  rc.tbl <- suppressMessages(read_delim(readcounts.file.ZOTU, delim = "\t")) %>%
    mutate(SampleID = as.character(SampleID))
  readcount.mat <- rc.tbl %>%
    select(-1) %>%
    as.matrix() %>%
    t() %>%
    set_colnames(rc.tbl$`#OTU ID`)
  metadata.tbl <- left_join(metadata.tbl,
                            data.frame(n_vsearch_reads_ZOTU = rowSums(readcount.mat),
                                       SampleID = rownames(readcount.mat)),
                            by = "SampleID")
  cat("Writing to file ", metadata.file, "\n")
  write_delim(metadata.tbl, delim = "\t", file = metadata.file)
}
