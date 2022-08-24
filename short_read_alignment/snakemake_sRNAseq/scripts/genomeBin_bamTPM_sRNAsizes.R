#!/usr/bin/env Rscript

# This R script is called by ../Snakefile
# Calculate per-window coverage normalised by total genome-wide coverage (TPM) for each chromosome
# for each sRNA size class                                              

# Usage: 
# ./genomeBin_bamTPM_sRNAsizes.R Col_0_sRNA_SRR1042171 t2t-col.20210610 both 10000"

args <- commandArgs(trailingOnly = TRUE)
sampleName <- args[1]
refbase <- args[2]
align <- args[3]
genomeBinSize <- as.integer(args[4])

if(floor(log10(genomeBinSize)) + 1 < 4) {
  genomeBinName <- paste0(genomeBinSize, "bp")
} else if(floor(log10(genomeBinSize)) + 1 >= 4 &
          floor(log10(genomeBinSize)) + 1 <= 6) {
  genomeBinName <- paste0(genomeBinSize/1e3, "kb")
} else if(floor(log10(genomeBinSize)) + 1 >= 7) {
  genomeBinName <- paste0(genomeBinSize/1e6, "Mb")
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

alnGR <- GRanges(readGAlignments(sampleFile))

## Calculate per-million scaling factor and windowed TPM
## based on alignments of sRNAs of all size classes,
## See https://rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
chr_RPK <- lapply(seq_along(chrs), function(x) {
  # Define adjacent windows
  binSeq <- seq(from = 1, to = chrLens[x], by = genomeBinSize)
  binIR <- IRanges(start = binSeq,
                   width = genomeBinSize)
  binIR <- binIR[-length(binIR)]
  binIR <- append(binIR,
                  IRanges(start = binSeq[length(binSeq)],
                          end = chrLens[x]))
  binGR <- GRanges(seqnames = chrs[x],
                   ranges = binIR,
                   strand = "*")

  # Count reads overlapping each bin
  bin_reads <- countOverlaps(query = binGR,
                             subject = alnGR,
                             type = "any",
                             ignore.strand = TRUE)
  # Divide read counts by the length of each bin in kb (reads per kilobase, RPK)
  bin_reads / (width(binGR)/1e+03)
})

# Sum all the RPK values and divide this number by 1e6
# (the "per million" scaling factor)
sample_RPK <- do.call(sum, chr_RPK)
RPKPM_scaling_factor <- sample_RPK/1e+06

# For each sRNA size class, calcualte TPM (transcripts per kilobase per million) at each genomic coordinate
# Divide each RPK value by RPKPM_scaling_factor to give TPM
TPM_sRNAsizes_df_list <- lapply(range(width(alnGR))[1]:range(width(alnGR))[2], function(w) {
  alnGR_sRNAsize <- alnGR[width(alnGR) == w]
  ## Calculate per-million scaling factor and windowed TPM
  ## based on alignments of sRNAs of all size classes,
  ## See https://rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
  chr_sRNAsize_RPK <- lapply(seq_along(chrs), function(x) {
    # Define adjacent windows
    binSeq <- seq(from = 1, to = chrLens[x], by = genomeBinSize)
    binIR <- IRanges(start = binSeq,
                     width = genomeBinSize)
    binIR <- binIR[-length(binIR)]
    binIR <- append(binIR,
                    IRanges(start = binSeq[length(binSeq)],
                            end = chrLens[x]))
    binGR <- GRanges(seqnames = chrs[x],
                     ranges = binIR,
                     strand = "*")
  
    # Count reads overlapping each bin
    bin_reads <- countOverlaps(query = binGR,
                               subject = alnGR_sRNAsize,
                               type = "any",
                               ignore.strand = TRUE)
    # Divide read counts by the length of each bin in kb (reads per kilobase, RPK)
    bin_reads / (width(binGR)/1e+03)
  })
  chr_sRNAsize_TPM <- lapply(seq_along(chr_sRNAsize_RPK), function(x) {
    chr_sRNAsize_RPK[[x]]/RPKPM_scaling_factor
  })
  chr_sRNAsize_TPM_df_list <- lapply(seq_along(chr_sRNAsize_TPM), function(x) {
    # Define adjacent windows
    binSeq <- seq(from = 1, to = chrLens[x], by = genomeBinSize)
    binIR <- IRanges(start = binSeq,
                     width = genomeBinSize)
    binIR <- binIR[-length(binIR)]
    binIR <- append(binIR,
                    IRanges(start = binSeq[length(binSeq)],
                            end = chrLens[x]))
    binGR <- GRanges(seqnames = chrs[x],
                     ranges = binIR,
                     strand = "*")
    data.frame(chr = chrs[x],
               window = start(binGR),
               cumwindow = start(binGR) + sumchr[x],
               cov = chr_sRNAsize_TPM[[x]])
  })
  return(do.call(rbind, chr_sRNAsize_TPM_df_list))
})
#}, mc.cores = length(range(width(alnGR))[1]:range(width(alnGR))[2]), mc.preschedule = F)

TPM_total <- do.call(sum, lapply(seq_along(TPM_sRNAsizes_df_list),
                            function(w) { TPM_sRNAsizes_df_list[[w]]$cov
                          }) )
print("Sum of TPM across all alignment sizes =")
print(TPM_total)

foreach(w = which(range(width(alnGR))[1]:range(width(alnGR))[2] %in% sRNAsizes)) %dopar% {
  print(paste0((range(width(alnGR))[1]:range(width(alnGR))[2])[w], "-nt sRNAs"))
  write.table(TPM_sRNAsizes_df_list[[w]],
              file = paste0("mapped/", align, "/tsv/",
                            sampleName, "_MappedOn_", refbase, "_", align, "_",
                            (range(width(alnGR))[1]:range(width(alnGR))[2])[w], "nt_sort_TPM_binSize", genomeBinName, ".tsv"),
              sep = "\t", quote = F, row.names = F, col.names = T)
}
