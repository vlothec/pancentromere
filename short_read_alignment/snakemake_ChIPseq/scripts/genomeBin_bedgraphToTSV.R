#!/home/ajt200/miniconda3/envs/R-4.0.0/bin/Rscript

# This R script is called by ../Snakefile
# Convert deepTools bamCoverage-generated *{genomeBinName}.bedgraph files into TSV files
# These TSV files can be imported into R for calculating and plotting log2(ChIP/control) chromosome-scale profiles 

# Usage:
# ./genomeBin_bedgraphToTSV.R WT_DMC1_V5_Rep1_ChIP TAIR10_chr_all both 10000

#sampleName <- "WT_DMC1_V5_Rep1_ChIP"
#refbase <- "TAIR10_chr_all"
#align <- "both"
#genomeBinSize <- 10000

args <- commandArgs(trailingOnly = T)
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
library(parallel)
library(plyr)
library(data.table)
#library(varhandle)
#library(zoo)

# Genomic definitions
fai <- read.table(paste0("data/index/", refbase, ".fa.fai"), header = F)
if(!grepl("Chr", fai[,1][1])) {
  chrs <- paste0("Chr", fai[,1])[1:5]
} else {
  chrs <- fai[,1][1:5]
}
chrLens <- fai[,2][1:5]

# Make chromosomal coordinates cumulative
# such that the first coordinate of Chr2 is
# equal to the last coordinate of Chr1 + 1
sumchr <- cumsum(c(0, chrLens))
print(sumchr)

# Load bedgraph
sampleProfile <- read.table(paste0("mapped/", align, "/bg/",
                                   sampleName, "_MappedOn_", refbase, "_lowXM_",
                                   align, "_sort_norm_binSize", genomeBinName, ".bedgraph"))
if(!grepl("Chr", fai[,1][1])) {
  sampleProfile$V1 <- paste0("Chr", sampleProfile$V1)
}
sampleProfile <- sampleProfile[sampleProfile$V1 %in% chrs,]
# Rows where the difference between end and start coordinates is > genomeBinSize
sampleProfile_bigWins <- sampleProfile[sampleProfile$V3-sampleProfile$V2 > genomeBinSize,]
# Rows where the difference between end and start coordinates is == genomeBinSize
sampleProfile <- sampleProfile[sampleProfile$V3-sampleProfile$V2 == genomeBinSize,]

# Create a list of big windows, each split into windows of genomeBinSize,
# or < genomeBinSize if at chromosome end
sampleProfile_bigWinsList <- mclapply(seq_along(1:dim(sampleProfile_bigWins)[1]), function(x) {
  bigWinsSplit <- seq(from = sampleProfile_bigWins[x,]$V2,
                      to = sampleProfile_bigWins[x,]$V3,
                      by = genomeBinSize)

  if(bigWinsSplit[length(bigWinsSplit)] < sampleProfile_bigWins[x,]$V3) {
    data.frame(V1 = as.character(sampleProfile_bigWins[x,]$V1),
               V2 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)],
                                 bigWinsSplit[length(bigWinsSplit)])),
               V3 = as.integer(c(bigWinsSplit[-length(bigWinsSplit)]+genomeBinSize,
                                 sampleProfile_bigWins[x,]$V3)),
               V4 = as.numeric(sampleProfile_bigWins[x,]$V4))
  } else if (bigWinsSplit[length(bigWinsSplit)] == sampleProfile_bigWins[x,]$V3) {
    data.frame(V1 = as.character(sampleProfile_bigWins[x,]$V1),
               V2 = as.integer(bigWinsSplit[-length(bigWinsSplit)]),
               V3 = as.integer(bigWinsSplit[-length(bigWinsSplit)]+genomeBinSize),
               V4 = as.numeric(sampleProfile_bigWins[x,]$V4))
  }
}, mc.cores = detectCores())

sampleProfile_bigWinsDT <- rbindlist(sampleProfile_bigWinsList)
sampleProfile <- rbind.fill(sampleProfile, sampleProfile_bigWinsDT)
sampleProfile <- sampleProfile[order(sampleProfile$V1, sampleProfile$V2),]

chrLenValsList <- mclapply(seq_along(chrs), function (x) {
  chrProfileSample <- sampleProfile[sampleProfile$V1 == chrs[x],]
  if(chrProfileSample[dim(chrProfileSample)[1],]$V3 < chrLens[x]) {
    data.frame(V1 = chrs[x],
               V2 = as.integer(chrProfileSample[dim(chrProfileSample)[1],]$V3),
               V3 = as.integer(chrLens[x]),
               V4 = as.numeric(chrProfileSample[dim(chrProfileSample)[1],]$V4))
  }
}, mc.cores = length(chrs))
sampleProfile_chrLenValsDT <- rbindlist(chrLenValsList)
sampleProfile <- rbind.fill(sampleProfile, sampleProfile_chrLenValsDT)
sampleProfile <- sampleProfile[order(sampleProfile$V1, sampleProfile$V2),]

chrCumWindowList <- lapply(seq_along(chrs), function(x) {
  chrProfileSample <- sampleProfile[sampleProfile$V1 == chrs[x],]
  chrProfileSample$V2+1+sumchr[x]
})

sampleProfile <- data.frame(chr = as.character(sampleProfile$V1),
                            window = as.integer(sampleProfile$V2+1),
                            cumwindow = as.integer(unlist(chrCumWindowList)),
                            cov = as.numeric(sampleProfile$V4))
write.table(sampleProfile,
            file = paste0("mapped/", align, "/tsv/",
                          sampleName, "_MappedOn_", refbase,
                          "_lowXM_", align, "_sort_norm_binSize", genomeBinName, ".tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)
