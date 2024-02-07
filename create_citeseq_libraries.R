.libPaths("/ihome/acillo/arc85/Rlibs_Mar_2023")
library(tidyverse)

cli <- commandArgs(trailingOnly=TRUE)
args <- strsplit(cli,"=",fixed=T)

## Read in samples
dat <- readr::read_csv("samples.csv")

## Create data for libraries.csv for input to cellranger count
dat <- dat %>%
  filter(sample==args[[1]])

to_write <- dat %>%
  select(fastq_citeseq,sample) %>%
  mutate(library_type="Antibody Capture")

colnames(to_write)[1] <- "fastqs"

readr::write_csv(to_write,file=paste("citeseq_libraries/",
  to_write$sample,
  "_CITEseq_library.csv",sep=""))