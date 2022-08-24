#!/usr/bin/env Rscript

# Compare average CEN180 metrics (HOR membership and divergence) in regions flanking
# centromeric ATHILA and matched centromeric random loci

# Usage:
# conda activate R-4.1.2
# ./66Atha_CENATHILA_flanks_CEN180_metrics_v260522_permTestAllList_hpc.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 1000 1e4 Flanks
# conda deactivate

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
flankSize <- as.numeric(args[2])
perms <- as.numeric(args[3])
regionName <- args[4]


# Set minimum possible P-value for permutation test result with
# perms sets of random loci
minPval <- 1 - ( (perms - 1) / perms)

if(floor(log10(flankSize)) + 1 < 4) {
  flankName <- paste0(flankSize, "bp")
} else if(floor(log10(flankSize)) + 1 >= 4 &
          floor(log10(flankSize)) + 1 <= 6) {
  flankName <- paste0(flankSize/1e3, "kb")
} else if(floor(log10(flankSize)) + 1 >= 7) {
  flankName <- paste0(flankSize/1e6, "Mb")
}
flankNamePlot <- paste0(c(strsplit(flankName, split = "")[[1]][1:(length(strsplit(flankName, split = "")[[1]])-2)],
                          "-",
                          strsplit(flankName, split = "")[[1]][(length(strsplit(flankName, split = "")[[1]])-1):(length(strsplit(flankName, split = "")[[1]]))]),
                          collapse = "")

options(stringsAsFactors = F)
source("/home/ajt200/rds/hpc-work/pancentromere/annotation/ATHILA/66Atha_CENATHILA_flanks_CEN180_metrics_v260522_permTestAllList_hpc_function.R")
suppressMessages(library(data.table, quietly = T))
suppressMessages(library(GenomicRanges, quietly = T))
suppressMessages(library(segmentSeq, quietly = T))
suppressMessages(library(parallel, quietly = T))
suppressMessages(library(dplyr, quietly = T))
suppressMessages(library(scales, quietly = T))
suppressMessages(library(pals, quietly = T))
suppressMessages(library(RColorBrewer, quietly = T))
suppressMessages(library(ggplot2, quietly = T))
suppressMessages(library(cowplot, quietly = T))
suppressMessages(library(ggbeeswarm, quietly = T))
suppressMessages(library(viridis, quietly = T))
suppressMessages(library(doParallel, quietly = T))
suppressMessages(library(doRNG, quietly = T))
suppressMessages(library(doFuture, quietly = T))
suppressMessages(library(snow, quietly = T))
suppressMessages(library(Rmpi, quietly = T))
suppressMessages(library(doMPI, quietly = T))
suppressMessages(library(iterators, quietly = T))
suppressMessages(library(plotrix, quietly = T))

# Create and register a doParallel parallel backend
registerDoParallel()
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())
options(future.globals.maxSize = +Inf)

# Define chunkSize so that each cluster worker gets a single "task chunk"
# (i.e. one task chunk covering multiple loop iterations (rows of binDF)),
# which is faster than if each cluster worker gets multiple task chunks
# (i.e., one task chunk per loop iteration (row of binDF))
chunkSize <- ceiling(perms / getDoParWorkers())
#initWorkers <- function() options(scipen = 999, stringsAsFactors = F)
mpiopts <- list(chunkSize = chunkSize)

plotDir <- paste0("ATHILA/plots/")
plotDirAllMetrics <- paste0(plotDir, "AllMetrics/")
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))
system(paste0("[ -d ", plotDirAllMetrics, " ] || mkdir -p ", plotDirAllMetrics))

# Load and define accession names
acc_full <- system("ls /home/ajt200/rds/hpc-work/pancentromere/annotation/CEN180/repeats/*.fa*", intern = T)
acc_full <- gsub("/home/ajt200/rds/hpc-work/pancentromere/annotation/CEN180/repeats/cen180.consensus.repetitiveness", "", acc_full)
acc_full <- acc_full[-grep("SUPER_", acc_full)]
acc_uniq <- unique( gsub("\\.fa\\..+", "", acc_full) )
acc_uniq_len <- NULL
for(i in acc_uniq) {
  acc_uniq_len <- c(acc_uniq_len, length( grep(i , acc_full)) )
}
acc <- acc_uniq[which(acc_uniq_len == length(chrName))]

chrs_list <- foreach(x = 1:length(acc), .inorder = T) %dopar% {
  acc_chrs <- read.table(paste0("/home/ajt200/rds/hpc-work/pancentromere/assemblies/",
                                acc[x], ".fa.fai"),
                         header = F)[,1]
  acc_chrs <- gsub("_RagTag_RagTag", "", acc_chrs)
  acc_chrs <- gsub("chr", "Chr", acc_chrs)
  acc_chrs <- gsub("SUPER_", "Chr", acc_chrs)
  acc_chrs <- acc_chrs[which(acc_chrs %in% chrName)]
  acc_chrs[sort.int(acc_chrs, index.return = T)$ix]
}

chrLens_list <- foreach(x = 1:length(acc), .inorder = T) %dopar% {
  acc_chrs <- read.table(paste0("/home/ajt200/rds/hpc-work/pancentromere/assemblies/",
                                acc[x], ".fa.fai"),
                         header = F)[,1]
  acc_chrs <- gsub("_RagTag_RagTag", "", acc_chrs)
  acc_chrs <- gsub("chr", "Chr", acc_chrs)
  acc_chrs <- gsub("SUPER_", "Chr", acc_chrs)
  acc_chrLens <- read.table(paste0("/home/ajt200/rds/hpc-work/pancentromere/assemblies/",
                                   acc[x], ".fa.fai"),
                            header = F)[,2]
  acc_chrLens <- acc_chrLens[which(acc_chrs %in% chrName)]
  acc_chrs <- acc_chrs[which(acc_chrs %in% chrName)]
  print(acc_chrs)
  print(acc_chrLens)
  acc_chrLens <- acc_chrLens[sort.int(acc_chrs, index.return = T)$ix]
  acc_chrs <- acc_chrs[sort.int(acc_chrs, index.return = T)$ix]
  print(acc_chrs)
  print(acc_chrLens)
  acc_chrLens
}

# Load CEN coordinates for each accession
CEN <- read.csv(paste0("/home/ajt200/rds/hpc-work/pancentromere/centromeric_coordinates/",
                       "centromere_manual_EDTA4_fa.csv"),
                header = T)
CEN$fasta.name <- gsub(".fa", "", CEN$fasta.name)
CEN_GR_list <- foreach(x = 1:length(acc), .inorder = T) %dopar% {
  acc_CEN <- CEN[grep(acc[x], CEN$fasta.name),]
  GRanges(seqnames = acc_CEN$chr,
          ranges = IRanges(start = acc_CEN$start,
                           end = acc_CEN$end),
          strand = "*",
          acc = acc_CEN$fasta.name,
          region = acc_CEN$region)
}

# Load CEN180 for each accession
CEN180_list <- foreach(x = 1:length(acc), .inorder = T) %dopar% {
  tab_list <- lapply(1:length(chrName), function(y) {
    read.csv(paste0("/home/ajt200/rds/hpc-work/pancentromere/annotation/CEN180/repeats/",
                    "cen180.consensus.repetitiveness", acc[x], ".fa.", chrName[y], ".csv"),
             header = T)
  })
  if(length(chrName) > 1) {
    tab <- dplyr::bind_rows(tab_list)
  } else {
    tab <- tab_list[[1]]
  }
  tab <- tab[which(tab$class == "aTha178"),]
  colnames(tab)[10] <- "chr"
  tab$chr <- gsub("_RagTag_RagTag", "", tab$chr)
  tab$chr <- gsub("chr", "Chr", tab$chr)
  tab$chr <- gsub("SUPER_", "Chr", tab$chr)
  tab$assembly <- gsub("\\.fasta", "", tab$assembly)
  tab$assembly <- gsub("\\.fa", "", tab$assembly)
  tab <- tab[
             with( tab, order(chr, start, end) ),
            ]
  tab
}

CEN180_GR_list <- foreach(x = 1:length(acc), .inorder = T) %dopar% { 
  GRanges(seqnames = as.character(CEN180_list[[x]]$chr),
          ranges = IRanges(start = as.integer(CEN180_list[[x]]$start),
                           end = as.integer(CEN180_list[[x]]$end)),
          strand = as.character(CEN180_list[[x]]$strand),
          acc = as.character(CEN180_list[[x]]$assembly),
          HORlengthsSum = as.numeric(CEN180_list[[x]]$HORlengthsSum),
          HORcount = as.numeric(CEN180_list[[x]]$HORcount),
          WeightedConsensusScore = as.numeric(CEN180_list[[x]]$weighted.consensus.score),
          EditDistance = as.numeric(CEN180_list[[x]]$edit.distance))
}

# Load CENATHILA for each accession
CENATHILA_list <- foreach(x = 1:length(acc), .inorder = T) %dopar% {
  tab <- read.table(paste0("/home/ajt200/rds/hpc-work/pancentromere/annotation/ATHILA/ATHILA/",
                           acc[x], "/CENATHILA_in_", acc[x], "_",
                           paste0(chrName, collapse = "_"), ".bed"),
                    header = F)
  tab$V2 <- tab$V2+1
  colnames(tab) <- c("chr", "start", "end", "name", "phylo", "strand")
  tab <- tab[
             with( tab, order(chr, start, end) ),
            ]
  tab
}

CENATHILA_DF_phylo <- dplyr::bind_rows(CENATHILA_list) 
phylo <- sort(unique(CENATHILA_DF_phylo$phylo))

phylo_colFun <- cols25(n = 18)[-c(7:16)][1:length(phylo)]
stopifnot(length(phylo_colFun) == length(phylo))
names(phylo_colFun) <- phylo

CENATHILA_GR_list <- foreach(x = 1:length(acc), .inorder = T) %dopar% {
  GRanges(seqnames = as.character(CENATHILA_list[[x]]$chr),
          ranges = IRanges(start = as.integer(CENATHILA_list[[x]]$start),
                           end = as.integer(CENATHILA_list[[x]]$end)),
          strand = as.character(CENATHILA_list[[x]]$strand),
          phylo = as.character(CENATHILA_list[[x]]$phylo))
}


# Define function to select randomly positioned loci of the same
# width distribution as CENgapAllAthila_bed
ranLocStartSelect <- function(coordinates, n) {
  sample(x = coordinates,
         size = n,
         replace = FALSE)
}

# Disable scientific notation (e.g., 59000000 rather than 5.9e+07)
options(scipen = 100)

# Function to define, for each accession, "perms" sets of centromeric random loci (acc_CENranLoc_GR)
# of the same number and width distribution as acc_CENATHILA_GR
defineCENranLoc <- function(acc_idx, chrs_list, chrLens_list, CEN_GR_list, CENATHILA_GR_list, seed) {
  print(acc[acc_idx])

  acc_chrs <- chrs_list[[acc_idx]]
  acc_chrLens <- chrLens_list[[acc_idx]]
  acc_CEN_GR <- CEN_GR_list[[acc_idx]]
  acc_CENATHILA_GR <- CENATHILA_GR_list[[acc_idx]]

  # Apply ranLocStartSelect() on a per-chromosome basis so that
  # acc_CENranLoc_GR contains the same number of loci per chromosome as acc_CENATHILA_GR
  acc_CENranLoc_GR <- GRanges()
  for(j in 1:length(acc_chrs)) {
    print(acc_chrs[j])

    chr_acc_CEN_GR <- acc_CEN_GR[seqnames(acc_CEN_GR) == acc_chrs[j]]

    chr_acc_CENATHILA_GR <- acc_CENATHILA_GR[seqnames(acc_CENATHILA_GR) == acc_chrs[j]]
    if(length(chr_acc_CENATHILA_GR) > 0) {
      set.seed(seed + 1e6)
      chr_acc_CENranLoc_Start <- ranLocStartSelect(coordinates = unlist(lapply(seq_along(chr_acc_CEN_GR), function(x) {
                                                                          ( start(chr_acc_CEN_GR[x]) + max(width(chr_acc_CENATHILA_GR)) ) :
                                                                          ( end(chr_acc_CEN_GR[x]) - max(width(chr_acc_CENATHILA_GR)) )
                                                                        })),
                                                   n = length(chr_acc_CENATHILA_GR))
      chr_acc_CENranLoc_GR <- GRanges(seqnames = acc_chrs[j],
                                      ranges = IRanges(start = chr_acc_CENranLoc_Start,
                                                       width = width(chr_acc_CENATHILA_GR)),
                                      strand = strand(chr_acc_CENATHILA_GR),
                                      phylo = as.character(chr_acc_CENATHILA_GR$phylo))
      acc_CENranLoc_GR <- append(acc_CENranLoc_GR, chr_acc_CENranLoc_GR)
    }
  }
  stopifnot(identical(width(acc_CENranLoc_GR), width(acc_CENATHILA_GR)))
  acc_CENranLoc_GR
}

# Define seed so that random selections are reproducible
set.seed(76492749)
CENranLoc_GR_acc_list <- foreach(acc_idx = 1:length(acc), .inorder = T) %do% {
  acc_perms_GR_list <- foreach(seed = iter(1:perms),
#                               .options.mpi = mpiopts,
                               .multicombine = T,
                               .maxcombine = perms+1e1,
                               .inorder = F) %dorng% {
    defineCENranLoc(acc_idx = acc_idx,
                    chrs_list = chrs_list,
                    chrLens_list = chrLens_list,
                    CEN_GR_list = CEN_GR_list,
                    CENATHILA_GR_list = CENATHILA_GR_list,
                    seed = seed)
  }
  acc_perms_GR_list
}

##Adapted from promoters() in GenomicRanges to extract regions relative to TTS
TTSplus <- function (x, upstream = 100, downstream = 1000, ...)
{
    if (!isSingleNumber(upstream))
        stop("'upstream' must be a single integer")
    if (!is.integer(upstream))
        upstream <- as.numeric(upstream)
    if (!isSingleNumber(downstream))
        stop("'downstream' must be a single integer")
    if (!is.integer(downstream))
        downstream <- as.numeric(downstream)
    if (downstream < 0)
        stop("'downstream' must be an integer >= 0")
#    if (upstream < 0 | downstream < 0)
#        stop("'upstream' and 'downstream' must be integers >= 0")
    if (any(strand(x) == "*"))
        warning("'*' ranges were treated as '+'")
    on_plus <- which(strand(x) == "+" | strand(x) == "*")
    on_plus_TTS <- end(x)[on_plus]
    start(x)[on_plus] <- on_plus_TTS - upstream
    end(x)[on_plus] <- on_plus_TTS + downstream
    on_minus <- which(strand(x) == "-")
    on_minus_TTS <- start(x)[on_minus]
    end(x)[on_minus] <- on_minus_TTS + upstream
    start(x)[on_minus] <- on_minus_TTS - downstream
    x
}

# Function to calculate mean CEN180 metrics for regions
# upstream and downstream of CENATHILA and CENranLoc bodies
CEN180metricsAtCENATHILA <- function(acc_idx, CEN180_GR, CENATHILA_GR, featureName) {
  print(acc[acc_idx])

  # Get flankSize bp upstream of start coordinates
  CENATHILA_up_GR <- promoters(CENATHILA_GR, upstream = flankSize, downstream = 0)
  CENATHILA_down_GR <- TTSplus(CENATHILA_GR, upstream = -1, downstream = flankSize)
  stopifnot(identical(CENATHILA_GR$phylo, CENATHILA_up_GR$phylo))
  stopifnot(identical(CENATHILA_GR$phylo, CENATHILA_down_GR$phylo))
  stopifnot(identical(seqnames(CENATHILA_GR), seqnames(CENATHILA_up_GR)))
  stopifnot(identical(seqnames(CENATHILA_GR), seqnames(CENATHILA_down_GR)))
  stopifnot(identical(strand(CENATHILA_GR), strand(CENATHILA_up_GR)))
  stopifnot(identical(strand(CENATHILA_GR), strand(CENATHILA_down_GR)))

  CENATHILA_CEN180_metrics <- data.frame()
  for(i in 1:length(chrName)) {
    print(chrName[i])

    CEN180_GR_chr <- CEN180_GR[seqnames(CEN180_GR) == chrName[i]]
    CENATHILA_GR_chr <- CENATHILA_GR[seqnames(CENATHILA_GR) == chrName[i]]
    CENATHILA_up_GR_chr <- CENATHILA_up_GR[seqnames(CENATHILA_up_GR) == chrName[i]]
    CENATHILA_down_GR_chr <- CENATHILA_down_GR[seqnames(CENATHILA_down_GR) == chrName[i]]
    stopifnot(identical(CENATHILA_GR_chr$phylo, CENATHILA_up_GR_chr$phylo))
    stopifnot(identical(CENATHILA_GR_chr$phylo, CENATHILA_down_GR_chr$phylo))
    stopifnot(identical(strand(CENATHILA_GR_chr), strand(CENATHILA_up_GR_chr)))
    stopifnot(identical(strand(CENATHILA_GR_chr), strand(CENATHILA_down_GR_chr)))

    if(length(CENATHILA_GR_chr) > 0) {
      CENATHILA_up_CEN180 <- getOverlaps(coordinates = CENATHILA_up_GR_chr,
                                         segments = CEN180_GR_chr,
                                         overlapType = "overlapping",
                                         whichOverlaps = TRUE,
                                         ignoreStrand = TRUE)
      CENATHILA_up_CEN180_HORlengthsSum <- sapply(1:length(CENATHILA_up_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_up_CEN180[[x]]]$HORlengthsSum, na.rm = T)
      })
      CENATHILA_up_CEN180_HORcount <- sapply(1:length(CENATHILA_up_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_up_CEN180[[x]]]$HORcount, na.rm = T)
      })
      CENATHILA_up_CEN180_WeightedConsensusScore <- sapply(1:length(CENATHILA_up_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_up_CEN180[[x]]]$WeightedConsensusScore, na.rm = T)
      })
      CENATHILA_up_CEN180_EditDistance <- sapply(1:length(CENATHILA_up_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_up_CEN180[[x]]]$EditDistance, na.rm = T)
      })

      CENATHILA_down_CEN180 <- getOverlaps(coordinates = CENATHILA_down_GR_chr,
                                           segments = CEN180_GR_chr,
                                           overlapType = "overlapping",
                                           whichOverlaps = TRUE,
                                           ignoreStrand = TRUE)
      CENATHILA_down_CEN180_HORlengthsSum <- sapply(1:length(CENATHILA_down_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_down_CEN180[[x]]]$HORlengthsSum, na.rm = T)
      })
      CENATHILA_down_CEN180_HORcount <- sapply(1:length(CENATHILA_down_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_down_CEN180[[x]]]$HORcount, na.rm = T)
      })
      CENATHILA_down_CEN180_WeightedConsensusScore <- sapply(1:length(CENATHILA_down_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_down_CEN180[[x]]]$WeightedConsensusScore, na.rm = T)
      })
      CENATHILA_down_CEN180_EditDistance <- sapply(1:length(CENATHILA_down_CEN180), function(x) {
        mean(CEN180_GR_chr[CENATHILA_down_CEN180[[x]]]$EditDistance, na.rm = T)
      })

      # Don't concatenate CENATHILA_up_GR_chr and CENATHILA_down_GR_chr;
      # better to keep start and end coordinates and annotate with "Upstream" or "Downstream" instead
      ##CENATHILA_reg_GR_chr <- c(CENATHILA_up_GR_chr, CENATHILA_down_GR_chr)
      stopifnot(length(CENATHILA_up_GR_chr) == length(CENATHILA_down_GR_chr))
      CENATHILA_reg_GR_chr <- c(CENATHILA_GR_chr, CENATHILA_GR_chr)
      CENATHILA_chr <- data.frame(CENATHILA_reg_GR_chr,
                                  feature = rep(featureName, length(CENATHILA_reg_GR_chr)),
                                  accession = rep(acc[acc_idx], length(CENATHILA_reg_GR_chr)),
                                  region = rep(c("Upstream", "Downstream"), each = length(CENATHILA_GR_chr)),
                                  HORlengthsSum = c(CENATHILA_up_CEN180_HORlengthsSum, CENATHILA_down_CEN180_HORlengthsSum),
                                  HORcount = c(CENATHILA_up_CEN180_HORcount, CENATHILA_down_CEN180_HORcount),
                                  WeightedConsensusScore = c(CENATHILA_up_CEN180_WeightedConsensusScore, CENATHILA_down_CEN180_WeightedConsensusScore),
                                  EditDistance = c(CENATHILA_up_CEN180_EditDistance, CENATHILA_down_CEN180_EditDistance))
      colnames(CENATHILA_chr)[1] <- "chr"
      CENATHILA_chr$chr <- as.character(CENATHILA_chr$chr)
      CENATHILA_chr$strand <- as.character(CENATHILA_chr$strand)
      colnames(CENATHILA_chr)[which(colnames(CENATHILA_chr) == "phylo")] <- "Family"
      colnames(CENATHILA_chr)[which(colnames(CENATHILA_chr) == "feature")] <- "Feature"
      colnames(CENATHILA_chr)[which(colnames(CENATHILA_chr) == "region")] <- "Region"

      CENATHILA_CEN180_metrics <- rbind(CENATHILA_CEN180_metrics, CENATHILA_chr)
    } 
  }
  CENATHILA_CEN180_metrics
}

# For each acc, apply CEN180metricsAtCENATHILA() function
# CENATHILA
CENATHILA_CEN180_metrics_list <- foreach(acc_idx = 1:length(acc), .inorder = T) %dopar% {
  CEN180metricsAtCENATHILA(acc_idx = acc_idx,
                           CEN180_GR = CEN180_GR_list[[acc_idx]],
                           CENATHILA_GR = CENATHILA_GR_list[[acc_idx]],
                           featureName = "CENATHILA")
}

# CENranLoc
CENranLoc_CEN180_metrics_list <- foreach(acc_idx = 1:length(acc), .inorder = T) %do% {
  foreach(x = iter(1:perms),
#          .options.mpi = mpiopts,
          .multicombine = T,
          .maxcombine = perms+1e1,
          .inorder = F) %dopar% {
    CEN180metricsAtCENATHILA(acc_idx = acc_idx,
                             CEN180_GR = CEN180_GR_list[[acc_idx]],
                             CENATHILA_GR = CENranLoc_GR_acc_list[[acc_idx]][[x]],
                             featureName = "CENranLoc")
  }
}


# HORlengthsSum
permTestAllList_HORlengthsSum <- permTestAllList(
                                                 CENATHILA_CEN180_metrics_list = CENATHILA_CEN180_metrics_list,
                                                 CENranLoc_CEN180_metrics_list = CENranLoc_CEN180_metrics_list,
                                                 region_name = regionName,
                                                 metric_name = "HORlengthsSum"
                                                )
permTestAllList_HORlengthsSum_permDistDF <- data.frame(
                                                       Accession = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                         rep(permTestAllList_HORlengthsSum[[y]]@accession,
                                                             times = length(permTestAllList_HORlengthsSum[[y]]@permDist))
                                                       })),
                                                       Family = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                         rep(permTestAllList_HORlengthsSum[[y]]@fam,
                                                             times = length(permTestAllList_HORlengthsSum[[y]]@permDist))
                                                       })),
                                                       Metric = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                         rep(permTestAllList_HORlengthsSum[[y]]@metric,
                                                             times = length(permTestAllList_HORlengthsSum[[y]]@permDist))
                                                       })),
                                                       Region = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                         rep(permTestAllList_HORlengthsSum[[y]]@region,
                                                             times = length(permTestAllList_HORlengthsSum[[y]]@permDist))
                                                       })),
                                                       Permuted = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                         permTestAllList_HORlengthsSum[[y]]@permDist
                                                       }))
                                                      ) 
permTestAllList_HORlengthsSum_permDF <- data.frame(
                                                   Accession = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@accession
                                                   })),
                                                   Family = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@fam
                                                   })),
                                                   Metric = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@metric
                                                   })),
                                                   Region = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@region
                                                   })),
                                                   Observed = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@observed
                                                   })),
                                                   Expected = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@expected
                                                   })),
                                                   AlphaThreshold0.05 = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@alphaThreshold
                                                   })),
                                                   Pvalue = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@pval
                                                   })),
                                                   AlternativeHypothesis = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@alternative
                                                   })),
                                                   Features = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@features
                                                   })),
                                                   Log2ObservedExpected = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@log2obsexp
                                                   })),
                                                   Log2AlphaThreshold0.05 = unlist(lapply(1:length(permTestAllList_HORlengthsSum), function(y) {
                                                     permTestAllList_HORlengthsSum[[y]]@log2alpha
                                                   }))
                                                  ) 

# EditDistance
permTestAllList_EditDistance <- permTestAllList(
                                                CENATHILA_CEN180_metrics_list = CENATHILA_CEN180_metrics_list,
                                                CENranLoc_CEN180_metrics_list = CENranLoc_CEN180_metrics_list,
                                                region_name = regionName,
                                                metric_name = "EditDistance"
                                               )
permTestAllList_EditDistance_permDistDF <- data.frame(
                                                      Accession = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                        rep(permTestAllList_EditDistance[[y]]@accession,
                                                            times = length(permTestAllList_EditDistance[[y]]@permDist))
                                                      })),
                                                      Family = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                        rep(permTestAllList_EditDistance[[y]]@fam,
                                                            times = length(permTestAllList_EditDistance[[y]]@permDist))
                                                      })),
                                                      Metric = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                        rep(permTestAllList_EditDistance[[y]]@metric,
                                                            times = length(permTestAllList_EditDistance[[y]]@permDist))
                                                      })),
                                                      Region = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                        rep(permTestAllList_EditDistance[[y]]@region,
                                                            times = length(permTestAllList_EditDistance[[y]]@permDist))
                                                      })),
                                                      Permuted = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                        permTestAllList_EditDistance[[y]]@permDist
                                                      }))
                                                     ) 
permTestAllList_EditDistance_permDF <- data.frame(
                                                  Accession = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@accession
                                                  })),
                                                  Family = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@fam
                                                  })),
                                                  Metric = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@metric
                                                  })),
                                                  Region = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@region
                                                  })),
                                                  Observed = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@observed
                                                  })),
                                                  Expected = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@expected
                                                  })),
                                                  AlphaThreshold0.05 = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@alphaThreshold
                                                  })),
                                                  Pvalue = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@pval
                                                  })),
                                                  AlternativeHypothesis = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@alternative
                                                  })),
                                                  Features = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@features
                                                  })),
                                                  Log2ObservedExpected = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@log2obsexp
                                                  })),
                                                  Log2AlphaThreshold0.05 = unlist(lapply(1:length(permTestAllList_EditDistance), function(y) {
                                                    permTestAllList_EditDistance[[y]]@log2alpha
                                                  }))
                                                  ) 

## HORcount
#permTestAllList_HORcount <- permTestAllList(
#                                            CENATHILA_CEN180_metrics_list = CENATHILA_CEN180_metrics_list,
#                                            CENranLoc_CEN180_metrics_list = CENranLoc_CEN180_metrics_list,
#                                            region_name = regionName,
#                                            metric_name = "HORcount"
#                                           )
#permTestAllList_HORcount_permDistDF <- data.frame(
#                                                  Accession = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                    rep(permTestAllList_HORcount[[y]]@accession,
#                                                        times = length(permTestAllList_HORcount[[y]]@permDist))
#                                                  })),
#                                                  Family = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                    rep(permTestAllList_HORcount[[y]]@fam,
#                                                        times = length(permTestAllList_HORcount[[y]]@permDist))
#                                                  })),
#                                                  Metric = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                    rep(permTestAllList_HORcount[[y]]@metric,
#                                                        times = length(permTestAllList_HORcount[[y]]@permDist))
#                                                  })),
#                                                  Region = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                    rep(permTestAllList_HORcount[[y]]@region,
#                                                        times = length(permTestAllList_HORcount[[y]]@permDist))
#                                                  })),
#                                                  Permuted = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                    permTestAllList_HORcount[[y]]@permDist
#                                                  }))
#                                                 ) 
#permTestAllList_HORcount_permDF <- data.frame(
#                                              Accession = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@accession
#                                              })),
#                                              Family = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@fam
#                                              })),
#                                              Metric = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@metric
#                                              })),
#                                              Region = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@region
#                                              })),
#                                              Observed = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@observed
#                                              })),
#                                              Expected = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@expected
#                                              })),
#                                              AlphaThreshold0.05 = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@alphaThreshold
#                                              })),
#                                              Pvalue = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@pval
#                                              })),
#                                              AlternativeHypothesis = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@alternative
#                                              })),
#                                              Features = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@features
#                                              })),
#                                              Log2ObservedExpected = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@log2obsexp
#                                              })),
#                                              Log2AlphaThreshold0.05 = unlist(lapply(1:length(permTestAllList_HORcount), function(y) {
#                                                permTestAllList_HORcount[[y]]@log2alpha
#                                              }))
#                                             ) 
#
## WeightedConsensusScore
#permTestAllList_WeightedConsensusScore <- permTestAllList(
#                                                          CENATHILA_CEN180_metrics_list = CENATHILA_CEN180_metrics_list,
#                                                          CENranLoc_CEN180_metrics_list = CENranLoc_CEN180_metrics_list,
#                                                          region_name = regionName,
#                                                          metric_name = "WeightedConsensusScore"
#                                                         )
#permTestAllList_WeightedConsensusScore_permDistDF <- data.frame(
#                                                                Accession = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                                  rep(permTestAllList_WeightedConsensusScore[[y]]@accession,
#                                                                      times = length(permTestAllList_WeightedConsensusScore[[y]]@permDist))
#                                                                })),
#                                                                Family = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                                  rep(permTestAllList_WeightedConsensusScore[[y]]@fam,
#                                                                      times = length(permTestAllList_WeightedConsensusScore[[y]]@permDist))
#                                                                })),
#                                                                Metric = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                                  rep(permTestAllList_WeightedConsensusScore[[y]]@metric,
#                                                                      times = length(permTestAllList_WeightedConsensusScore[[y]]@permDist))
#                                                                })),
#                                                                Region = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                                  rep(permTestAllList_WeightedConsensusScore[[y]]@region,
#                                                                      times = length(permTestAllList_WeightedConsensusScore[[y]]@permDist))
#                                                                })),
#                                                                Permuted = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                                  permTestAllList_WeightedConsensusScore[[y]]@permDist
#                                                                }))
#                                                               ) 
#permTestAllList_WeightedConsensusScore_permDF <- data.frame(
#                                                            Accession = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@accession
#                                                            })),
#                                                            Family = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@fam
#                                                            })),
#                                                            Metric = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@metric
#                                                            })),
#                                                            Region = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@region
#                                                            })),
#                                                            Observed = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@observed
#                                                            })),
#                                                            Expected = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@expected
#                                                            })),
#                                                            AlphaThreshold0.05 = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@alphaThreshold
#                                                            })),
#                                                            Pvalue = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@pval
#                                                            })),
#                                                            AlternativeHypothesis = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@alternative
#                                                            })),
#                                                            Features = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@features
#                                                            })),
#                                                            Log2ObservedExpected = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@log2obsexp
#                                                            })),
#                                                            Log2AlphaThreshold0.05 = unlist(lapply(1:length(permTestAllList_WeightedConsensusScore), function(y) {
#                                                              permTestAllList_WeightedConsensusScore[[y]]@log2alpha
#                                                            }))
#                                                           ) 



save(permTestAllList_HORlengthsSum,
     file = paste0(plotDirAllMetrics,
                   "CENATHILA_", flankName, "_regions_CEN180_HORlengthsSum_all_accessions_combined_",
                   paste0(chrName, collapse = "_"), "_",
                   perms, "perms_",
                   regionName,
                   "_permObject.RData"))
load(file = paste0(plotDirAllMetrics,
                   "CENATHILA_", flankName, "_regions_CEN180_HORlengthsSum_all_accessions_combined_",
                   paste0(chrName, collapse = "_"), "_",
                   perms, "perms_",
                   regionName,
                   "_permObject.RData"))

save(permTestAllList_EditDistance,
     file = paste0(plotDirAllMetrics,
                   "CENATHILA_", flankName, "_regions_CEN180_EditDistance_all_accessions_combined_",
                   paste0(chrName, collapse = "_"), "_",
                   perms, "perms_",
                   regionName,
                   "_permObject.RData"))
load(file = paste0(plotDirAllMetrics,
                   "CENATHILA_", flankName, "_regions_CEN180_EditDistance_all_accessions_combined_",
                   paste0(chrName, collapse = "_"), "_",
                   perms, "perms_",
                   regionName,
                   "_permObject.RData"))

#save(permTestAllList_HORcount,
#     file = paste0(plotDirAllMetrics,
#                   "CENATHILA_", flankName, "_regions_CEN180_HORcount_all_accessions_combined_",
#                   paste0(chrName, collapse = "_"), "_",
#                   perms, "perms_",
#                   regionName,
#                   "_permObject.RData"))
#load(file = paste0(plotDirAllMetrics,
#                   "CENATHILA_", flankName, "_regions_CEN180_HORcount_all_accessions_combined_",
#                   paste0(chrName, collapse = "_"), "_",
#                   perms, "perms_",
#                   regionName,
#                   "_permObject.RData"))
#
#save(permTestAllList_WeightedConsensusScore,
#     file = paste0(plotDirAllMetrics,
#                   "CENATHILA_", flankName, "_regions_CEN180_WeightedConsensusScore_all_accessions_combined_",
#                   paste0(chrName, collapse = "_"), "_",
#                   perms, "perms_",
#                   regionName,
#                   "_permObject.RData"))
#load(file = paste0(plotDirAllMetrics,
#                   "CENATHILA_", flankName, "_regions_CEN180_WeightedConsensusScore_all_accessions_combined_",
#                   paste0(chrName, collapse = "_"), "_",
#                   perms, "perms_",
#                   regionName,
#                   "_permObject.RData"))


# Combine into one data.frame for plotting with ggplot2
combined_permDistDF <- rbind(
                             permTestAllList_HORlengthsSum_permDistDF,
                             permTestAllList_EditDistance_permDistDF
#                             permTestAllList_HORcount_permDistDF,
#                             permTestAllList_WeightedConsensusScore_permDistDF
                            )
combined_permDF <- rbind(
                         permTestAllList_HORlengthsSum_permDF,
                         permTestAllList_EditDistance_permDF
#                         permTestAllList_HORcount_permDF,
#                         permTestAllList_WeightedConsensusScore_permDF
                        )

write.table(combined_permDistDF,
            file = paste0(plotDirAllMetrics,
                          "CENATHILA_", flankName, "_regions_CEN180_AllMetrics_all_accessions_combined_",
                          paste0(chrName, collapse = "_"), "_",
                          perms, "perms_",
                          regionName,
                          "_permDistDF.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
combined_permDistDF <- fread(file = paste0(plotDirAllMetrics,
                                           "CENATHILA_", flankName, "_regions_CEN180_AllMetrics_all_accessions_combined_",
                                           paste0(chrName, collapse = "_"), "_",
                                           perms, "perms_",
                                           regionName,
                                           "_permDistDF.tsv"),
                             data.table = F)
write.table(combined_permDF,
            file = paste0(plotDirAllMetrics,
                          "CENATHILA_", flankName, "_regions_CEN180_AllMetrics_all_accessions_combined_",
                          paste0(chrName, collapse = "_"), "_",
                          perms, "perms_",
                          regionName,
                          "_permDF.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
combined_permDF <- fread(file = paste0(plotDirAllMetrics,
                                       "CENATHILA_", flankName, "_regions_CEN180_AllMetrics_all_accessions_combined_",
                                       paste0(chrName, collapse = "_"), "_",
                                       perms, "perms_",
                                       regionName,
                                       "_permDF.tsv"),
                         data.table = F)

combined_permDistDF$Accession <- factor(combined_permDistDF$Accession,
                                        levels = sort(unique(combined_permDistDF$Accession)))
combined_permDistDF$Family <- factor(combined_permDistDF$Family,
                                     levels = rev(sort(unique(combined_permDistDF$Family))))
combined_permDistDF$Metric <- factor(combined_permDistDF$Metric,
                                     levels =  unique(combined_permDistDF$Metric))
combined_permDistDF$Region <- factor(combined_permDistDF$Region,
                                     levels = unique(combined_permDistDF$Region))
combined_permDF$Accession <- factor(combined_permDF$Accession,
                                    levels = sort(unique(combined_permDF$Accession)))
combined_permDF$Family <- factor(combined_permDF$Family,
                                 levels = rev(sort(unique(combined_permDF$Family))))
combined_permDF$Metric <- factor(combined_permDF$Metric,
                                 levels =  unique(combined_permDF$Metric))
combined_permDF$Region <- factor(combined_permDF$Region,
                                 levels = unique(combined_permDF$Region))


# Define function to make colours transparent,
# to aid visibility where points overlap
makeTransparent <- function(thisColour, alpha = 230)
{
  newColour <- col2rgb(thisColour)
  apply(newColour, 2, function(x) {
    rgb(red = x[1], green = x[2], blue = x[3],
        alpha = alpha, maxColorValue = 255)
  })
}


vp_all <- ggplot(data = combined_permDistDF,
                 mapping = aes(x = Family,
                               y = Permuted)) +
  xlab("Family") +
  ylab(bquote("Mean" ~ italic("CEN178") ~ "metric in" ~ .(flankNamePlot) ~ "flanking regions")) +
  geom_violin(trim = F,
              scale = "count",
              colour = "grey70",
              fill = "grey70") +
  geom_point(data = combined_permDF,
             mapping = aes(x = Family,
                           y = AlphaThreshold0.05),
             position = position_nudge(x = -0.04),
             shape = 124, colour = makeTransparent("darkorange1"), size = 12) +
  geom_point(data = combined_permDF,
             mapping = aes(x = Family,
                           y = Observed),
             position = position_nudge(x = -0.04),
             shape = 124, colour = makeTransparent("dodgerblue2"), size = 12) +
  geom_point(data = combined_permDF,
             mapping = aes(x = Family,
                           y = Expected),
             position = position_nudge(x = -0.04),
             shape = 124, colour = makeTransparent("black"), size = 12) +

  coord_flip() +
  theme_bw() +
  theme(
#        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.ticks.y = element_line(size = 0.5, colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black", face = "italic", hjust = 0.5, vjust = 0.5, angle = 0),
        axis.title.y = element_text(size = 20, colour = "black"),
#        axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.ticks.x = element_line(size = 0.5, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black", hjust = 0.5, vjust = 0.5, angle = 0),
        axis.title.x = element_text(size = 20, colour = "black", margin = margin(t = 20)),
        strip.text.x = element_text(size = 20, colour = "black"),
        strip.text.y = element_text(size = 20, colour = "black"),
        strip.background = element_rect(fill = "white", colour = "white"),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        plot.title = element_text(size = 15, colour = "black", hjust = 0.5)) +
  ggtitle(bquote("Permutations =" ~
                 .(prettyNum(perms,
                             big.mark = ",",
                             trim = T)) ~
                 "sets of randomly positioned centromeric loci"))

vp_all <- vp_all +
  facet_grid(cols = vars(Metric), scales = "free_x")
ggsave(paste0(plotDirAllMetrics,
              "CENATHILA_", flankName, "_regions_CEN180_AllMetrics_all_accessions_violin_",
              paste0(chrName, collapse = "_"), "_",
              perms, "perms_",
              paste0(unique(as.vector(combined_permDF$Region)), collapse = "_"),
              ".pdf"),
       plot = vp_all,
       width = 4*length(unique(combined_permDF$Metric)), height = 8, limitsize = F)


# Plot histogram summaries

# HORlengthsSum
pdf(paste0(plotDirAllMetrics,
           "CENATHILA_", flankName, "_regions_CEN180_HORlengthsSum_ATHILA_all_accessions_histogram_",
           paste0(chrName, collapse = "_"), "_",
           perms, "perms_",
           paste0(unique(as.vector(combined_permDF$Region)), collapse = "_"),
           ".pdf"),
    height = 4.5, width = 5.0)
par(mar = c(3.1, 3.1, 4.1, 1.1),
    mgp = c(1.85, 0.75, 0))
# Disable scientific notation (e.g., 0.0001 rather than 1e-04)
options(scipen = 100)
# Calculate max density
maxDensityPlus <- max(density(permTestAllList_HORlengthsSum[[1]]@permDist)$y)*1.1
if(permTestAllList_HORlengthsSum[[1]]@alternative == "MoreThanRandom") {
  xlim <- c(pmin(min(permTestAllList_HORlengthsSum[[1]]@permDist)/1.1),
            pmax(permTestAllList_HORlengthsSum[[1]]@observed*1.1, permTestAllList_HORlengthsSum[[1]]@alphaThreshold*1.1))
  textX1 <- quantile(c(pretty(xlim)[1],
                       pretty(xlim)[length(pretty(xlim))]), 0.10)[[1]]
#  textX1 <- min(permTestAllList_HORlengthsSum[[1]]@permDist)/1.15
} else {
  xlim <- c(pmin(permTestAllList_HORlengthsSum[[1]]@observed/1.1),
            max(permTestAllList_HORlengthsSum[[1]]@permDist)*1.1)
  textX1 <- quantile(xlim, 0.75)[[1]]
#  textX1 <- min(permTestAllList_HORlengthsSum[[1]]@permDist)/1.15
}
hist(permTestAllList_HORlengthsSum[[1]]@permDist,
     breaks = 50,
     freq = FALSE,
     col = "grey70",
     border = NA,
     lwd = 2,
     xlim = c(pretty(xlim)[1],
              pretty(xlim)[length(pretty(xlim))]),
     ylim = c(0,
              maxDensityPlus),
     xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", main = "",
     axes = FALSE)
axis(side = 2,
     at = pretty(density(permTestAllList_HORlengthsSum[[1]]@permDist)$y),
     lwd = 1)
mtext(side = 2,
      text = "Density",
      line = 1.85)
axis(side = 1,
     at = pretty(xlim),
     lwd = 1)
mtext(side = 1,
      text = bquote("Mean" ~ italic("CEN178") ~ "HORlengthsSum in" ~ .(flankNamePlot) ~ "flanking regions"),
      line = 1.85)
titleText <- list(bquote("Centromeric" ~ italic(.(permTestAllList_HORlengthsSum[[1]]@fam))),
                  bquote(italic("P")*" = "*
                         .(as.character(permTestAllList_HORlengthsSum[[1]]@pval))),
#                         .(as.character(round(permTestAllList_HORlengthsSum[[1]]@pval,
#                                              digits = 6)))),
                  bquote("Permutations = "*.(prettyNum(length(permTestAllList_HORlengthsSum[[1]]@permDist),
                                                       big.mark = ",",
                                                       trim = T)) ~
                         "sets of randomly positioned centromeric loci"))
mtext(do.call(expression, titleText), side = 3, line = 3:1, cex = c(1, 0.7, 0.7))
lines(density(permTestAllList_HORlengthsSum[[1]]@permDist),
      col = "grey70",
      lwd = 1.5)
ablineclip(v = permTestAllList_HORlengthsSum[[1]]@expected,
           y1 = 0, y2 = maxDensityPlus*.92, lwd = 2.5)
ablineclip(v = permTestAllList_HORlengthsSum[[1]]@observed,
           y1 = 0, y2 = maxDensityPlus*.92, lwd = 2.5, col = "dodgerblue2")
ablineclip(v = permTestAllList_HORlengthsSum[[1]]@alphaThreshold,
           y1 = 0, y2 = maxDensityPlus*.92, lwd = 2.5, lty = 5, col = "darkorange1")
text(x = c(textX1,
           permTestAllList_HORlengthsSum[[1]]@expected,
           permTestAllList_HORlengthsSum[[1]]@observed,
           permTestAllList_HORlengthsSum[[1]]@alphaThreshold),
     y = c(maxDensityPlus*.95,
           maxDensityPlus,
           maxDensityPlus,
           maxDensityPlus*.95),
     labels = c("Permuted",
                "Expected",
                "Observed",
                expression(alpha*" = 0.05")),
     col = c("grey70",
             "black",
             "dodgerblue2",
             "darkorange1"),
     cex = 0.8)
dev.off()


# EditDistance
pdf(paste0(plotDirAllMetrics,
           "CENATHILA_", flankName, "_regions_CEN180_EditDistance_ATHILA_all_accessions_histogram_",
           paste0(chrName, collapse = "_"), "_",
           perms, "perms_",
           paste0(unique(as.vector(combined_permDF$Region)), collapse = "_"),
           ".pdf"),
    height = 4.5, width = 5.0)
par(mar = c(3.1, 3.1, 4.1, 1.1),
    mgp = c(1.85, 0.75, 0))
# Disable scientific notation (e.g., 0.0001 rather than 1e-04)
options(scipen = 100)
# Calculate max density
maxDensityPlus <- max(density(permTestAllList_EditDistance[[1]]@permDist)$y)*1.01
if(permTestAllList_EditDistance[[1]]@alternative == "MoreThanRandom") {
  xlim <- c(pmin(min(permTestAllList_EditDistance[[1]]@permDist)/1.01),
            pmax(permTestAllList_EditDistance[[1]]@observed*1.01, permTestAllList_EditDistance[[1]]@alphaThreshold*1.01))
  textX1 <- quantile(c(pretty(xlim)[1],
                       pretty(xlim)[length(pretty(xlim))]), 0.10)[[1]]
#  textX1 <- min(permTestAllList_EditDistance[[1]]@permDist)/1.15
} else {
  xlim <- c(pmin(permTestAllList_EditDistance[[1]]@observed/1.01),
            max(permTestAllList_EditDistance[[1]]@permDist)*1.01)
  textX1 <- quantile(xlim, 0.75)[[1]]
#  textX1 <- min(permTestAllList_EditDistance[[1]]@permDist)/1.15
}
hist(permTestAllList_EditDistance[[1]]@permDist,
     breaks = 50,
     freq = FALSE,
     col = "grey70",
     border = NA,
     lwd = 2,
     xlim = c(pretty(xlim)[1],
              pretty(xlim)[length(pretty(xlim))]),
     ylim = c(0,
              maxDensityPlus),
     xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", main = "",
     axes = FALSE)
axis(side = 2,
     at = pretty(density(permTestAllList_EditDistance[[1]]@permDist)$y),
     lwd = 1)
mtext(side = 2,
      text = "Density",
      line = 1.85)
axis(side = 1,
     at = pretty(xlim),
     lwd = 1)
mtext(side = 1,
      text = bquote("Mean" ~ italic("CEN178") ~ "EditDistance in" ~ .(flankNamePlot) ~ "flanking regions"),
      line = 1.85)
titleText <- list(bquote("Centromeric" ~ italic(.(permTestAllList_EditDistance[[1]]@fam))),
                  bquote(italic("P")*" = "*
                         .(as.character(permTestAllList_EditDistance[[1]]@pval))),
#                         .(as.character(round(permTestAllList_EditDistance[[1]]@pval,
#                                              digits = 6)))),
                  bquote("Permutations = "*.(prettyNum(length(permTestAllList_EditDistance[[1]]@permDist),
                                                       big.mark = ",",
                                                       trim = T)) ~
                         "sets of randomly positioned centromeric loci"))
mtext(do.call(expression, titleText), side = 3, line = 3:1, cex = c(1, 0.7, 0.7))
lines(density(permTestAllList_EditDistance[[1]]@permDist),
      col = "grey70",
      lwd = 1.5)
ablineclip(v = permTestAllList_EditDistance[[1]]@expected,
           y1 = 0, y2 = maxDensityPlus*.92, lwd = 2.5)
ablineclip(v = permTestAllList_EditDistance[[1]]@observed,
           y1 = 0, y2 = maxDensityPlus*.92, lwd = 2.5, col = "dodgerblue2")
ablineclip(v = permTestAllList_EditDistance[[1]]@alphaThreshold,
           y1 = 0, y2 = maxDensityPlus*.92, lwd = 2.5, lty = 5, col = "darkorange1")
text(x = c(textX1,
           permTestAllList_EditDistance[[1]]@expected,
           permTestAllList_EditDistance[[1]]@observed,
           permTestAllList_EditDistance[[1]]@alphaThreshold),
     y = c(maxDensityPlus*.95,
           maxDensityPlus,
           maxDensityPlus,
           maxDensityPlus*.95),
     labels = c("Permuted",
                "Expected",
                "Observed",
                expression(alpha*" = 0.05")),
     col = c("grey70",
             "black",
             "dodgerblue2",
             "darkorange1"),
     cex = 0.8)
dev.off()



print("warnings 1")
print(warnings())

## Shutdown the cluster and quit
##stopCluster(cl) # use if cl made with makeCluster() (i.e., using doFuture package)
#closeCluster(cl) # use if cl made with startMPIcluster() (i.e., using doMPI package - faster than doFuture)
#mpi.quit()
#
#print("warnings 2")
#print(warnings())

