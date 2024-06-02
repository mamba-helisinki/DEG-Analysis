library(DESeq2)
library(pheatmap)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(biomaRt)


countData <- read.table("C:/Users/HP/Desktop/GSE153104_all_counts.tsv/GSE153104_all_counts.xlsx", header=TRUE, row.names=1)

