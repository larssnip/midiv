#' @name demultiplex
#' @title De-multiplexing fastq files
#'
#' @description De-multiplexing Illumina data based on the extra forward barcode
#' used by MiDiv lab.
#'
#' @param sample.tbl Table with data for each sample, see Details below.
#' @param in.folder	Name of folder where raw fastq files are located.
#' @param out.folder Name of folder to output de-multiplexed fastq files.
#' @param compress.out Logical to indicate compressed output or not.
#' @param pattern The pattern to recognize the raw fastq files from other files
#'
#' @details The input \code{sample.tbl} must be a table (tibble or data.frame)
#' with one row for each sample. It must follow the MiDiv sample_table standard
#' format. The columns used by this function are:
#' * SampleID
#' * Rawfile_R1
#' * Rawfile_R2
#' * Barcode
#'
#' The SampleID should be a short text to uniquely identify each sample. It will
#' be used to form the new fastq filenames, so avoid using symbols not recommended
#' in filenames (e.g. space, slash). It is also wise to have a letter, not an integer,
#' as the first symbol in such SampleID texts.
#'
#' Demultiplexing means extracting subsets of reads from raw fastq files, and saving
#' each subset on new fastq files. Each subset is identified by a barcode
#' sequence, and this must be listed in the Barcode column. The Barcode
#' sequence is matched at the start of the R1-reads only.
#'
#' The original fastq files must all be in the \code{in.folder}. These files may
#' be compressed (.gz). The \code{pattern} is matched to all file names in this
#' folder, and should give a match against all raw fastq files. All new fastq
#' files written to \code{out.folder} will
#' also be compressed (.gz) if \code{compress.out = TRUE}.
#'
#' @return The function will output the de-multiplexed fastq-files to the
#' \code{out.folder}. The name of each file consists of the corresponding
#' SampleID text, with the \code{_R1.fastq}/\code{_R2.fastq} extension. If
#' \code{compress.out = TRUE} the extension \code{.gz} is also added.
#'
#' The function will return in R a table with the number of read-pairs for each
#' sample. You may then add this as a new column to the existing
#' \code{sample.tbl} by
#' \code{full_join(sample.tbl, demultiplex.tbl, by = "SampleID")}, where
#' \code{demultiplex.tbl} indicates the output from this function.
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
demultiplex <- function(sample.tbl, in.folder, out.folder, compress.out = FALSE,
                        pattern = "_R[12].fastq"){
  cat("De-multiplexing: ")
  in.folder <- normalizePath(in.folder)
  out.folder <- normalizePath(out.folder)
  out.ext <- c("_R1.fastq", "_R2.fastq")
  if(compress.out) out.ext <- str_c(out.ext, ".gz")
  cnames <- c("SampleID", "Rawfile_R1", "Rawfile_R2", "Barcode")
  if(sum(is.na(match(cnames, colnames(sample.tbl)))) > 0)
    stop("The sample.tbl must have columns named SampleID, Rawfile_R1, Rawfile_R2 and Barcode")
  utbl <- sample.tbl %>%
    select(Rawfile_R1, Rawfile_R2) %>%
    distinct()
  readpairs.tbl <- tibble(SampleID = sample.tbl$SampleID,
                          n_readpairs = 0)
  for(i in 1:nrow(utbl)) {
    cat("   Reading raw file", utbl$Rawfile_R1[i], "...")
    R1.tbl <- readFastq(file.path(in.folder, utbl$Rawfile_R1[i])) %>%
      rename(R1.Header = Header, R1.Sequence = Sequence, R1.Quality = Quality)
    cat("   Reading raw file", utbl$Rawfile_R2[i], "...")
    R2.tbl <- readFastq(file.path(in.folder, utbl$Rawfile_R2[i])) %>%
      rename(R2.Header = Header, R2.Sequence = Sequence, R2.Quality = Quality)
    tbl <- bind_cols(R1.tbl, R2.tbl)
    idx <- which(sample.tbl$Rawfile_R1 == utbl$Rawfile_R1[i] & sample.tbl$Rawfile_R2 == utbl$Rawfile_R2[i])
    cat(" found", length(idx), "samples with reads from these raw files:\n")
    for(j in 1:length(idx)){
      cat("      Barcode", sample.tbl$Barcode[idx[j]], "...\n")
      nc <- str_length(sample.tbl$Barcode[idx[j]])
      tbl0 <- tbl %>%
        filter(str_detect(R1.Sequence, str_c("^", sample.tbl$Barcode[idx[j]]))) %>%
        mutate(R1.Sequence = str_sub(R1.Sequence, start = nc + 1, end = -1)) %>%
        mutate(R1.Quality = str_sub(R1.Quality, start = nc + 1, end = -1))
      tbl0 %>% select(starts_with("R1")) %>%
        rename(Header = R1.Header, Sequence = R1.Sequence, Quality = R1.Quality) %>%
        writeFastq(out.file = file.path(out.folder, str_c(sample.tbl$SampleID[idx[j]], out.ext[1])))
      tbl0 %>% select(starts_with("R2")) %>%
        rename(Header = R2.Header, Sequence = R2.Sequence, Quality = R2.Quality) %>%
        writeFastq(out.file = file.path(out.folder, str_c(sample.tbl$SampleID[idx[j]], out.ext[2])))
      cat("         found", nrow(tbl0), "read-pairs\n")
      readpairs.tbl$n_readpairs[idx[j]] <- nrow(tbl0)
    }
  }
  return(readpairs.tbl)
}
