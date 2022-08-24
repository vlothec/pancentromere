#!/usr/bin/env Rscript

# This R script is called by ../Snakefile
# Calculate per-base coverage normalised by total genome-wide coverage (TPM) for each chromosome
# for all sRNAs and for each sRNA size class                                              

# Usage: 
# ./bin_bamTPM_sRNAsizes.R Col_0_sRNA_SRR1042171 t2t-col.20210610 both 1"

args <- commandArgs(trailingOnly = TRUE)
sampleName <- args[1]
refbase <- args[2]
align <- args[3]
binSize <- as.integer(args[4])

if(floor(log10(binSize)) + 1 < 4) {
  binName <- paste0(binSize, "bp")
} else if(floor(log10(binSize)) + 1 >= 4 &
          floor(log10(binSize)) + 1 <= 6) {
  binName <- paste0(binSize/1e3, "kb")
} else if(floor(log10(binSize)) + 1 >= 7) {
  binName <- paste0(binSize/1e6, "Mb")
}

options(stringsAsFactors = F)
library(yaml)
library(Rsamtools)
library(GenomicAlignments)
library(parallel)
library(doParallel)
registerDoParallel(cores = detectCores())
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

config <- read_yaml("config.yaml")
sRNAsizes <- config$sRNA_SIZES
#sRNAsizes <- c(18:26)

sampleFile <- paste0("mapped/", align, "/",
                     sampleName, "_MappedOn_", refbase, "_", align, "_sort.bam")

# Genomic definitions
fai <- read.table(paste0("data/index/", refbase, ".fa.fai"), header = F)
ignoreForNormalization <- unlist(strsplit(config$COVERAGE$ignoreForNormalization,
                                          split = " "))
fai <- fai[!(fai$V1 %in% ignoreForNormalization),]
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1])
} else {
  chrs <- fai[,1]
}
chrLens <- fai[,2]

# Make chromosomal coordinates cumulative
# such that the first coordinate of Chr2 is
# equal to the last coordinate of Chr1 + 1
sumchr <- cumsum(c(0, chrLens))
print(sumchr)

## Calculate TPM
## See https://rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
# Divide read counts by the length of each bin in kb (reads per kilobase, RPK)
rGA <- readGAlignments(sampleFile)
rGA_cov <- coverage(rGA)
chr_RPK <- lapply(chrs, function(x) {
  as.numeric(rGA_cov[names(rGA_cov) == x][[1]]) / (binSize/1e+03)
})
# Sum all the RPK values and divide this number by 1e6
# (the "per million" scaling factor)
sample_RPK <- do.call(sum, chr_RPK)
RPKPM_scaling_factor <- sample_RPK/1e+06

# For each sRNA size class, calcualte TPM (transcripts per kilobase per million) at each genomic coordinate
# Divide each RPK value by RPKPM_scaling_factor to give TPM
TPM_sRNAsizes_df_list <- lapply(range(width(rGA))[1]:range(width(rGA))[2], function(w) {
  rGA_sRNAsize <- rGA[width(rGA) == w]
  rGA_sRNAsize_cov <- coverage(rGA_sRNAsize)
  chr_sRNAsize_RPK <- lapply(chrs, function(x) {
    as.numeric(rGA_sRNAsize_cov[names(rGA_sRNAsize_cov) == x][[1]]) / (binSize/1e+03)
  })
  chr_sRNAsize_TPM <- lapply(seq_along(chr_sRNAsize_RPK), function(x) {
    chr_sRNAsize_RPK[[x]]/RPKPM_scaling_factor
  })
  chr_sRNAsize_TPM_df_list <- lapply(seq_along(chr_sRNAsize_TPM), function(x) {
    data.frame(chr = chrs[x],
               start = 0:(chrLens[x]-1),
               end = 1:chrLens[x],
               TPM = chr_sRNAsize_TPM[[x]])
  })
  return(do.call(rbind, chr_sRNAsize_TPM_df_list))
})
#}, mc.cores = length(range(width(rGA))[1]:range(width(rGA))[2]), mc.preschedule = F)

TPM_total <- do.call(sum, lapply(seq_along(TPM_sRNAsizes_df_list),
                            function(w) { TPM_sRNAsizes_df_list[[w]]$TPM
                          }) )
print("Sum of TPM across all alignment sizes =")
print(TPM_total)

foreach(w = which(range(width(rGA))[1]:range(width(rGA))[2] %in% sRNAsizes)) %dopar% {
#for(w in which(range(width(rGA))[1]:range(width(rGA))[2] %in% sRNAsizes)) {
  print(paste0((range(width(rGA))[1]:range(width(rGA))[2])[w], "-nt sRNAs"))
  write.table(TPM_sRNAsizes_df_list[[w]],
              file = paste0("mapped/", align, "/bg/",
                            sampleName, "_MappedOn_", refbase, "_", align, "_",
                            (range(width(rGA))[1]:range(width(rGA))[2])[w], "nt_sort_TPM.bedgraph"),
              sep = "\t", quote = F, row.names = F, col.names = F)
}
