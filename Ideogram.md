---
title: "Chromosome Ideogram"
author: "Elena Polonsky"
date: "2026-08-03"
output: word_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

Load required packages
```{r}
require(RIdeogram)
library(rtracklayer)
library(Biostrings)
library(GenomicRanges)
```

Package versions:

* RIdeogram v 0.2.2  
* rtracklayer v 1.62.0  
* Biostrings v 2.70.3  
* GenomicRanges v 1.54.1  

```{r}
# Load annotation and genome file
liftoff_annotation <- import("~/Desktop/consensus_lifted_annotation.gff3")
genome <- readDNAStringSet("~/Desktop/consensus.fasta")
```

```{r}

# Filter only three chromosomes to be in the ideogram and exclude the unplaced scaffold
chromosomes <- c(
    "OZ004311.1",
    "OZ004312.1",
    "OZ004313.1"
)

# Subset the genome to retain only the selected chromosomes
genome <- genome[chromosomes]

#  Extract only gene features from the annotation file
genes <- liftoff_annotation[liftoff_annotation$type == "gene"]

# Keep only genes located on the selected chromosomes.
genes <- genes[seqnames(genes) %in% chromosomes]

# Calculate the length of each chromosome in base pairs
chr_lengths <- width(genome)

# Create the karyotype data frame required by RIdeogram
karyotype <- data.frame(
    Chr = names(genome),
    Start = 0,
    End = chr_lengths
)
```

```{r}
# Define the genomic window size used to calculate gene density
# 1 Mb intervals are used

window <- 1000000

# Calculate gene density across each chromosome
# Each chromosome is divided into consecutive windows and the number of genes overlapping each window is counted

gene_density <- do.call(rbind, lapply(seq_along(genome), function(i){
  
    chr <- names(genome)[i]
    chr_length <- width(genome)[i]
    
    windows <- seq(1, chr_length, by = window)

    data.frame(
        Chr = chr,
        Start = windows,
        End = pmin(windows + window - 1, chr_length),
        Value = sapply(windows, function(start){
            end <- min(start + window - 1, chr_length)

            sum(
                seqnames(genes) == chr &
                start(genes) <= end &
                end(genes) >= start
            )
        })
    )
}))
```

```{r}
# Generate the chromosome ideogram using the chromosome coordinates
# Overlay the calculated gene density

ideogram(
    karyotype = karyotype,
    overlaid = gene_density
)

convertSVG("chromosome.svg", device = "png")
knitr::include_graphics("chromosome.png")
```
