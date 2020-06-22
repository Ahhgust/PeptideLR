#!/usr/bin/Rscript
# Written by August Woerner
# This takes in multiple files from peptideLR
# and combines the values to generate a naive likelihood ratio
# estimate...
# Taking the max LR may not be what you want to do, so be mindful!

suppressPackageStartupMessages(library(tidyverse))



argv <- commandArgs(trailingOnly=TRUE)

if (length(argv)==0) {
    stop("I need 1 or more files from peptideLR from the *same analysis* to continue")
}


read1file <- function(f) {
    
    tib <- readr::read_tsv(f, col_types=cols(), progress=FALSE)
    tib$Filename <- f
    return(tib)
}

# read all the files in the argv
combined <- lapply(argv, read1file) %>%
    dplyr::bind_rows()

if (! "experimentType"  %in% colnames(combined)) {
    stop("Unexpected input...")
}
    
# sanity check.
exps <- unique(combined$experimentType)
if (n_distinct(combined$experimentType) > 1) {
    write(unique(combined$experimentType), stderr())
    stop("Multiple experiment-types are present...\nThat's not good!")
}

allcols <- colnames(combined)
# column names may vary; these remain the same.
groups <- allcols[ allcols != "Denom" & allcols !="Num" & allcols != "Filename" & allcols != "Theta"]

group_by_at(combined, groups) %>%
    summarize(NumProduct=prod(Num),
              DenomProduct=prod(Denom),
              LR=NumProduct/DenomProduct) %>%
    ungroup() %>% arrange(desc(LR)) -> summ

cat( readr::format_tsv(summ) )



