#' @name demultiplex
#' @title De-multiplexing fastq files
#'
#' @description De-multiplexing Illumina data based on the extra forward barcode
#' used by MiDiv lab.
#'
#' @param metadata.tbl Table with data for each sample, see Details below.
#' @param in.folder	Name of folder where the SequencingRunID folder with raw fastq files are located.
#' @param out.folder Name of folder to output de-multiplexed fastq files.
#' @param compress.out Logical to indicate compressed output or not.
#' @param pattern The pattern to recognize the raw fastq files from other files
#' @param trim.primers Logical indicating if PCR-primers should be trimmed from start of R1 and R2 reads.
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
#' * Forward_primer
#' * Reverse_primer
#'
#' The ProjectID, SequencingRunID and SampleID should all be a short text
#' (SampleID may be just an integer). The names of the de-multiplexed fastq files will
#' follow the format: ProjectID_SequencingRunID_SampleID_Rx.fastq.gz, where x is 1 or 2,
#' so avoid using symbols not recommended in file names (e.g. space, slash).
#'
#' Note that SequencingRunID must be the name of the folder in which the Rawfile_R1
#' and Rawfile_R2 fastq files are found. The \code{in.folder} argument is the path
#' to the SequencingRunID folder. Example: Raw files are in the folder
#' /mnt/rawdata/illumina/20250101_testrun/. Then \code{in.folder} is "/mnt/rawdata/illumina"
#' and \code{SequencingRunID} is "20250101_testrun".
#'
#' De-multiplexing means extracting subsets of reads from raw fastq files, those
#' named in columns Rawfile_R1 and Rawfile_R2 (if single-end reads, only Rawfile_R1).
#' The subset of read-pairs for each sample is identified by a barcode
#' sequence, and this must be listed in the Barcode column. The Barcode
#' sequence is matched at the start of the R1-reads only.
#'
#' If \code{trim.primers=TRUE} the start of the R1 sequence is trimmed by the
#' length of Forward_primer, and the start of the R2 read trimmed by the length
#' of Reverse_primer. NOTE: There is no primer-matching here. No reads are discarded,
#' only trimmed by primer lengths.
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
#' @importFrom stringr str_c str_length str_detect
#' @importFrom dplyr %>% distinct select filter mutate rename bind_cols
#' @importFrom microseq readFastq writeFastq
#'
#' @export demultiplex
#'
demultiplex <- function(metadata.tbl, in.folder, out.folder, trim.primers = TRUE){
  cat("De-multiplexing:\n")
  in.folder <- normalizePath(in.folder)
  out.folder <- normalizePath(out.folder)
  cnames <- c("ProjectID", "SequencingRunID", "SampleID", "Rawfile_R1", "Rawfile_R2", "Barcode", "Forward_primer", "Reverse_primer")
  if(sum(is.na(match(cnames, colnames(metadata.tbl)))) > 0)
    stop("The metadata.tbl must have columns named ProjectID, SequencingRunID, SampleID, Rawfile_R1, Rawfile_R2, Barcode, Forward_primer and Reverse_primer")
  utbl <- metadata.tbl %>%
    select(SequencingRunID, Rawfile_R1, Rawfile_R2) %>%
    distinct()
  readpairs.tbl <- tibble(ProjectID = metadata.tbl$ProjectID,
                          SequencingRunID = metadata.tbl$SequencingRunID,
                          SampleID = metadata.tbl$SampleID,
                          n_readpairs = 0)
  for(i in 1:nrow(utbl)) {
    cat("   Reading raw file", utbl$Rawfile_R1[i], "...\n")
    R1.tbl <- readFastq(file.path(in.folder, utbl$SequencingRunID[i], utbl$Rawfile_R1[i])) %>%
      rename(R1.Header = Header, R1.Sequence = Sequence, R1.Quality = Quality)
    cat("   Reading raw file", utbl$Rawfile_R2[i], "...\n")
    R2.tbl <- readFastq(file.path(in.folder, utbl$SequencingRunID[i], utbl$Rawfile_R2[i])) %>%
      rename(R2.Header = Header, R2.Sequence = Sequence, R2.Quality = Quality)
    tbl <- bind_cols(R1.tbl, R2.tbl)
    idx <- which(metadata.tbl$Rawfile_R1 == utbl$Rawfile_R1[i] & metadata.tbl$Rawfile_R2 == utbl$Rawfile_R2[i])
    cat("   found", length(idx), "samples with reads from these raw files:\n")
    for(j in 1:length(idx)){
      cat("      Barcode", metadata.tbl$Barcode[idx[j]], "...\n")
      nc <- str_length(metadata.tbl$Barcode[idx[j]])
      nf <- str_length(metadata.tbl$Forward_primer[idx[j]])
      nr <- str_length(metadata.tbl$Reverse_primer[idx[j]])
      R1.tbl <- R1.tbl |>
        mutate(start2 = str_sub(R1.Sequence, 1, nc + 2))
      M <- str_locate(R1.tbl$start2, metadata.tbl$Barcode[idx[j]])
      rr <- which(!is.na(M[,1]))
      tbl0 <- tbl %>%
        slice(rr) %>%
        mutate(R1.Sequence = str_sub(R1.Sequence, start = M[rr,2] + 1, end = -1)) %>%
        mutate(R1.Quality = str_sub(R1.Quality, start = M[rr,2] + 1, end = -1))
        # filter(str_detect(start2, metadata.tbl$Barcode[idx[j]]))) %>%
        # mutate(R1.Sequence = str_sub(R1.Sequence, start = nc + 1, end = -1)) %>%
        # mutate(R1.Quality = str_sub(R1.Quality, start = nc + 1, end = -1))
      if(trim.primers){
        tbl0 <- tbl0 %>%
          mutate(R1.Sequence = str_sub(R1.Sequence, start = nf + 1, end = -1)) %>%
          mutate(R1.Quality = str_sub(R1.Quality, start = nf + 1, end = -1)) %>%
          mutate(R2.Sequence = str_sub(R2.Sequence, start = nr + 1, end = -1)) %>%
          mutate(R2.Quality = str_sub(R2.Quality, start = nr + 1, end = -1))
      }
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
