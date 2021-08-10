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

# lazy man's flag processing
minContamination<- 0.00
if (argv[[1]] == '-c') {
    if (length(argv) < 3) {
        stop("Usage: -c dropin file1 file2 ... You specific -c and nothing else?!")
    }
    
    minContamination <- as.numeric(argv[[2]])
    if (minContamination > 1.0) {
        stop("Usage: -c dropin where dropin is <1")
    }
    # remove flag
    argv <- argv[-(1:2)]
}


read1file <- function(f) {
    
    tib <- readr::read_tsv(f, col_types=cols(), progress=FALSE, guess_max=10000)
    if (nrow(tib)>0) {
        tib$Filename <- f
    } else {
        add_column(tib, "Filename")
    }
    return(tib)
}

# read all the files in the argv
combined <- lapply(argv, read1file) %>%
    dplyr::bind_rows()



if (! "experimentType"  %in% colnames(combined)) {
    stop("Unexpected input...")
}
    

if (n_distinct(combined$versionNumber) > 1) {
    write(unique(combined$versionNumber), stderr())
    stop("Multiple version numbers are present...\nThat's not good!")
}

if (n_distinct(combined$experimentType) > 1) {
    write(unique(combined$experimentType), stderr())
    stop("Multiple experiment types...\nThat's not good!")
}

allcols <- colnames(combined)


if ("RMP_theta" %in% allcols) {
    
    # RMP processing.
    combined %>%
        filter(Dropin >= minContamination) %>%
        group_by(
            Dropin,
            Population,
            versionNumber,
            randomSeed,
            LikelihoodTransform
        ) %>%
        summarize(
            Npeps=sum(Npeps),
            RMP_theta=prod(RMP_theta),
            RMP_theta_0.025=prod(RMP_theta_0.025),
            RMP_theta_0.975=prod(RMP_theta_0.975),
            RMP_naive=prod(RMP_naive),
            RMP_naive_0.025=prod(RMP_naive_0.025),
            RMP_naive_0.975=prod(RMP_naive_0.975)
            ) %>%
        ungroup()  -> summ

} else {

    if (combined$experimentType[[1]] == "LR_Known_AllSingleSource") {
        combined$NTotal<-1 # definitional. number of contribs is 1.
                                        # higher mixture degrees have this as a variable...
        allcols <- colnames(combined)
    }
    

  # LR processing...
  # column names may vary; these remain the same.
    groups <- allcols[ allcols !="Likelihood" & allcols != "Filename" & allcols != "Theta"]
    
  # deprecated in tidyverse (group_by_at)... which is just stupid. 
# get likelihoods at the level of the genome,
# for every combination of drop-in/drop-out
    group_by_at(combined, groups) %>%
        filter(Dropin_Rate >= minContamination) %>%
        summarize(Likelihood=prod(Likelihood),
                  NCombined=n(),
                  ) %>%
        ungroup() %>% arrange(desc(Likelihood)) -> likes
    
    
   if (n_distinct(likes$NCombined) > 1) {
        stop("Heterogenous combinations detected. Some of your LR files appear to be truncated..?!")    }

    groups <- groups[ groups!='Dropin_Rate' & groups != 'Dropout_Rate']
    # extract out the MLE over drop-in and drop-out
    maxlikes <- group_by_at(likes, groups) %>%
        top_n(1, Likelihood) %>%
        ungroup()

    # TODO: flag in making plots...
    if (FALSE) {
        filter(likes,ID=="Random") %>%
            ggplot(aes(x=Dropin_Rate, y=log10(Likelihood), color=factor(NTotal))) +
            geom_point() +
            facet_grid(~Dropout_Rate)
        
        
        ord <- unique(maxlikes$ID)
        maxlikes$ID <- factor(maxlikes$ID, levels=ord)
        
        filter(maxlikes, NKnown < NTotal | NTotal==1) %>%
            ggplot(aes(x=ID, y=Likelihood,
                       color=factor(Dropin_Rate)
                       )) +
            geom_point() +
            scale_y_log10()+
            coord_flip() +
            facet_wrap(~NTotal, scales='free') +
            theme_bw(base_size=20) +
            labs(y="Likelihood", x="ID", color="Drop-in\nrate") -> pl

        ggsave("Likelihoods.png", pl, height=8, width=15)
    }
    
    rando <- filter(maxlikes, ID=="Random" | ID=="RandomPerson")
    # combined to give likelihood ratios
    # (only wrt to a denominator of ALL unknowns...)
    left_join(
        maxlikes,
        select(rando, NullLike=Likelihood, NTotal, Population),
        by=c("NTotal", "Population")) -> summ
        
    if (FALSE) {
        filter(summ, NKnown < NTotal | NTotal==1) %>%
            ggplot(aes(x=ID, y=log10(Likelihood/NullLike), color=factor(Dropin_Rate))) +
            geom_point() +
            coord_flip() +
            facet_wrap(~NTotal, scales='free') +
            theme_bw(base_size=20) +
            labs(y="Log10 likelihood ratio", x="ID", color="Drop-in\nrate") -> pl


        ggsave("LRsVsUnknowns.png", pl, height=8, width=15)
        
        filter(summ, Likelihood>=NullLike) %>%
            ggplot(aes(x=ID, y=log10(Likelihood/NullLike), color=factor(Dropin_Rate))) +
            geom_point() +
            coord_flip() +
            facet_wrap(~NTotal, scales='free') +
            theme_bw(base_size=20) +
            labs(y="Log10 likelihood ratio", x="ID", color="Drop-in\nrate") -> pl

        ggsave("AllPositiveLRs.png", pl, height=12, width=15)
    }
    
}

cat( readr::format_tsv(summ) )



