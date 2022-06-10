#' @name demultiplex
#' @title De-multiplexing fastq files
#'
#' @description De-multiplexing Illumina data based on the extra forward barcode
#' used by MiDiv lab.
#'
#' @param metadata.tbl Table with data for each sample, see Details below.
#' @param in.folder	Name of folder where raw fastq files are located.
#' @param out.folder Name of folder to output de-multiplexed fastq files.
#' @param compress.out Logical to indicate compressed output or not.
#' @param pattern The pattern to recognize the raw fastq files from other files
#'
#' @details The input \code{metadata.tbl} must be a table (tibble or data.frame)
#' with one row for each sample. It must follow the MiDiv metadata table standard
#' format. The columns used by this function are:
#' * ProjectID
#' * SequencingRunID
#' * SampleID
#' * Rawfile_R1
#' * Rawfile_R2
#' * Barcode
#'
#' The ProjectID, SequencingRunID and SampleID should all be a short text
#' (sampleID may be just an integer). The names of the de-multiplexed fastq files will
#' follow the format: ProjectID_SequencingRunID_SampleID_Rx.fastq.gz, where x is 1 or 2,
#' so avoid using symbols not recommended in filenames (e.g. space, slash).
#'
#' De-multiplexing means extracting subsets of reads from raw fastq files, those
#' named in columns Rawfile_R1 and Rawfile_R2 (if single-end reads, only Rawfile_R1).
#' The subset of read-pairs for each sample is identified by a barcode
#' sequence, and this must be listed in the Barcode column. The Barcode
#' sequence is matched at the start of the R1-reads only.
#'
#' The files listed in Rawfile_R1 and Rawfile_R2 must all be in the \code{in.folder}.
#' These files may be compressed (.gz).
#'
#' @return The function will output the de-multiplexed fastq-files to the
#' \code{out.folder}. The name of each file consists of the corresponding
#' ProjectID_SequencingRunID_SampleID, with the extensions \code{_R1.fastq.gz}
#' or \code{_R2.fastq.gz}.
#'
#' The function will return in R a table with the number of read-pairs for each
#' sample. You may then add this as a new column to the existing
#' \code{metadata.tbl} by
#' \code{full_join(metadata.tbl, demultiplex.tbl, by = c("ProjectID", "SequencingRunID", "SampleID")},
#' where \code{demultiplex.tbl} indicates the output from this function.
#'
#'
#' @author Lars Snipen.
#'
#' @importFrom stringr str_c str_length
#' @importFrom dplyr filter %>% if_else select mutate rename bind_cols
#' @importFrom microseq readFastq writeFastq
#'
#' @export readFasta
#' @export writeFasta
#'
demultiplex <- function(metadata.tbl, in.folder, out.folder){
  cat("De-multiplexing: ")
  in.folder <- normalizePath(in.folder)
  out.folder <- normalizePath(out.folder)
  cnames <- c("ProjectID", "SequencingRunID", "SampleID", "Rawfile_R1", "Rawfile_R2", "Barcode")
  if(sum(is.na(match(cnames, colnames(metadata.tbl)))) > 0)
    stop("The metadata.tbl must have columns named ProjectID, SequencingRunID, SampleID, Rawfile_R1, Rawfile_R2 and Barcode")
  utbl <- metadata.tbl %>%
    select(Rawfile_R1, Rawfile_R2) %>%
    distinct()
  readpairs.tbl <- tibble(ProjectID = metadata.tbl$ProjectID,
                          SequencingRunID = metadata.tbl$SequencingRunID,
                          SampleID = metadata.tbl$SampleID,
                          n_readpairs = 0)
  for(i in 1:nrow(utbl)) {
    cat("   Reading raw file", utbl$Rawfile_R1[i], "...")
    R1.tbl <- readFastq(file.path(in.folder, utbl$Rawfile_R1[i])) %>%
      rename(R1.Header = Header, R1.Sequence = Sequence, R1.Quality = Quality)
    cat("   Reading raw file", utbl$Rawfile_R2[i], "...")
    R2.tbl <- readFastq(file.path(in.folder, utbl$Rawfile_R2[i])) %>%
      rename(R2.Header = Header, R2.Sequence = Sequence, R2.Quality = Quality)
    tbl <- bind_cols(R1.tbl, R2.tbl)
    idx <- which(metadata.tbl$Rawfile_R1 == utbl$Rawfile_R1[i] & metadata.tbl$Rawfile_R2 == utbl$Rawfile_R2[i])
    cat(" found", length(idx), "samples with reads from these raw files:\n")
    for(j in 1:length(idx)){
      cat("      Barcode", metadata.tbl$Barcode[idx[j]], "...\n")
      nc <- str_length(metadata.tbl$Barcode[idx[j]])
      tbl0 <- tbl %>%
        filter(str_detect(R1.Sequence, str_c("^", metadata.tbl$Barcode[idx[j]]))) %>%
        mutate(R1.Sequence = str_sub(R1.Sequence, start = nc + 1, end = -1)) %>%
        mutate(R1.Quality = str_sub(R1.Quality, start = nc + 1, end = -1))
      tbl0 %>% select(starts_with("R1")) %>%
        rename(Header = R1.Header, Sequence = R1.Sequence, Quality = R1.Quality) %>%
        writeFastq(out.file = file.path(out.folder, str_c(metadata.tbl$ProjectID[idx[j]], "_",
                                                          metadata.tbl$SequencingRunID[idx[j]], "_",
                                                          metadata.tbl$SampleID[idx[j]], "_R1.fastq.gz")))
      tbl0 %>% select(starts_with("R2")) %>%
        rename(Header = R2.Header, Sequence = R2.Sequence, Quality = R2.Quality) %>%
        writeFastq(out.file = file.path(out.folder, str_c(metadata.tbl$ProjectID[idx[j]], "_",
                                                          metadata.tbl$SequencingRunID[idx[j]], "_",
                                                          metadata.tbl$SampleID[idx[j]], "_R2.fastq.gz")))
      cat("         found", nrow(tbl0), "read-pairs\n")
      readpairs.tbl$n_readpairs[idx[j]] <- nrow(tbl0)
    }
  }
  return(readpairs.tbl)
}
