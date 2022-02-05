#' @name demultiplex
#' @title Demultiplexing fastq files
#'
#' @description Demultiplexing Illumina data based on the extra forward barcode
#' used by MiDiv lab.
#'
#' @param sample.tbl Table with data for each sample, see Details below.
#' @param in.folder	Name of folder where raw fastq files are located.
#' @param out.folder Name of older to output demultiplexed fastq files.
#' @param fastq.tag Short text (file extension) to identify fastq files from other files.
#' @param R1.tag Short text to identify R1 (forward) reads files.
#' @param R2.tag Short text to identify R2 (reverse) reads files.
#' @param compress.out Logical to indicate compressed output or not.
#'
#' @details The input \code{sample.tbl} must be a table (tibble or data.frame)
#' with one row for each sample. It must contain the columns named (exactly)
#' SampleID, BarcodeSequence and FastqFileTag, but may have other columns as well.
#' Below the 3 required columns are explained.
#'
#' The SampleID should be a short text to uniquely identify each sample. It will
#' be used to form the new fastq filenames, so avoid using symbols not recommended
#' in filenames (e.g. space, slash). It is also wise to have a letter, not an integer,
#' as the first symbol in such SampleID texts.
#'
#' Demultiplexing means extracting subsets of reads from fastq files, and saving
#' each subset on new fastq files. Each subset is identified by a barcode
#' sequence, and this must be listed in the BarcodeSequence column. The Barcode
#' sequence is matched at the start of the R1-reads only.
#'
#' The original fastq files must all be in the \code{in.folder}. The \code{fastq.tag}
#' is used to separate these files from other files that may be present in \code{in.folder}.
#' The FastqFileTag must be a short text that uniquely identifies each raw fastaq
#' file in \code{in.folder}, i.e. it must be a part of the filename that is
#' unique to each file-pair (could be the entire filename).
#'
#' The fastq files in \code{in.folder} may be compressed (.gz). All new fastq
#' files written to out.folder will be compressed (.gz) if compress.out=TRUE.
#'
#' @return The function will output the demultiplexed fastq-files to the \code{out.folder}.
#' The name of each file consists of the corresponding SampleID text, followed
#' by the \code{R1.tag}/\code{R2.tag}, followed by the \code{fastq.tag}. If
#' \code{compress.out = TRUE} then the extension `.gz` is also added.
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
demultiplex <- function (sample.tbl, in.folder, out.folder, fastq.tag = "fastq",
                         R1.tag = "R1", R2.tag = "R2", compress.out = TRUE){
  cat("Demultiplexing:\n")
  out.ext <- if_else(compress.out, str_c(".", fastq.tag, ".gz"), str_c(".", fastq.tag))
  cnames <- c("SampleID", "BarcodeSequence", "FastqFileTag")
  if(sum(is.na(match(cnames, colnames(sample.tbl)))) > 0)
    stop("The sample.tbl must have columns named SampleID, BarcodeSequence and FastqFileTag")
  utag <- unique(sample.tbl$FastqFileTag)
  raw.files <- list.files(normalizePath(in.folder), pattern = fastq.tag, full.names = T)
  cat("found", length(raw.files), "fastq files in in.folder\n")
  if(length(raw.files) == 0)
    stop("No raw fastq files to demultiplex!")
  if(length(raw.files)%%2 != 0)
    stop("Must have pairs of fastq files!")
  for(i in 1:length(utag)) {
    cat("   Reading raw fastq files matching", utag[i], "...")
    raw <- raw.files[str_detect(raw.files, utag[i])]
    if(length(raw) == 0)
      stop("The FastqFileTag", utag[i], "has no matches in the fastq filenames")
    if(length(raw) == 1)
      stop("The FastqFileTag", utag[i], "produces only 1 match in the fastq filenames, should be 2")
    if(length(raw) > 2)
      stop("The FastqFileTag", utag[i], "produces more than 2 matches in the fastq filenames, should be 2")
    R1.raw <- raw[str_detect(raw, R1.tag)]
    if(length(R1.raw) != 1)
      stop("No match to", R1.tag, "in filenames", raw)
    R1.tbl <- readFastq(R1.raw) %>%
      rename(R1.Header = Header, R1.Sequence = Sequence, R1.Quality = Quality)
    R2.raw <- raw[str_detect(raw, R2.tag)]
    if(length(R2.raw) != 1)
      stop("No match to", R2.tag, "in filenames", raw)
    R2.tbl <- readFastq(R2.raw) %>%
      rename(R2.Header = Header, R2.Sequence = Sequence, R2.Quality = Quality)
    tbl <- bind_cols(R1.tbl, R2.tbl)
    idx <- which(str_detect(sample.tbl$FastqFileTag, utag[i]))
    cat(" found", length(idx), "samples:\n")
    for(j in 1:length(idx)){
      cat("      Barcode", sample.tbl$BarcodeSequence[idx[j]], "...")
      nc <- str_length(sample.tbl$BarcodeSequence[idx[j]])
      tbl0 <- tbl %>%
        filter(str_detect(R1.Sequence, str_c("^", sample.tbl$BarcodeSequence[idx[j]]))) %>%
        mutate(R1.Sequence = str_sub(R1.Sequence, start = nc + 1, end = -1)) %>%
        mutate(R1.Quality = str_sub(R1.Quality, start = nc + 1, end = -1))
      tbl0 %>% select(starts_with("R1")) %>%
        rename(Header = R1.Header, Sequence = R1.Sequence, Quality = R1.Quality) %>%
        writeFastq(out.file = file.path(normalizePath(out.folder), str_c(sample.tbl$SampleID[idx[j]], "_", R1.tag, out.ext)))
      tbl0 %>% select(starts_with("R2")) %>%
        rename(Header = R2.Header, Sequence = R2.Sequence, Quality = R2.Quality) %>%
        writeFastq(out.file = file.path(normalizePath(out.folder), str_c(sample.tbl$SampleID[idx[j]], "_", R2.tag, out.ext)))
      cat("found", nrow(tbl0), "read-pairs\n")
    }
  }
  return(NULL)
}
