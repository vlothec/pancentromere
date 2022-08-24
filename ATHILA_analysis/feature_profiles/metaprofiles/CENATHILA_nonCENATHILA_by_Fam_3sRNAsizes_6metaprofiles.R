#!/usr/bin/env Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 23.03.2022

# Calculate and plot metaprofiles of sRNA size classes
# around CENATHILA and nonCENATHIA within a given family

# Usage:
# ./CENATHILA_nonCENATHILA_by_Fam_3sRNAsizes_6metaprofiles.R 'Chr1,Chr2,Chr3,Chr4,Chr5' both 180 2000 2000 2000 2kb 10 10 10bp 10bp '0.02,0.96' 'AGO5_D7_sRNA_Rep1,Input_D7_sRNA_Rep1' 'miRNAseq_seedling_floral_Bradamante_Gutzat_2022_bioRxiv/snakemake_miRNAseq_,miRNAseq_seedling_floral_Bradamante_Gutzat_2022_bioRxiv/snakemake_miRNAseq_' 'AGO5 D7 Rep1,Input D7 Rep1' 'deepskyblue3,gold,seagreen3,orange1,mediumpurple4,orange3' '21,22,24' ATHILA1 t2t-col.20210610

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#align <- "both"
#bodyLength <- 180
#ATHILA_bodyLength <- 2000
#TEsf_bodyLength <- 2000
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#binSize <- 10
#ATHILA_binSize <- 10
#binName <- "10bp"
#ATHILA_binName <- "10bp"
## top left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.96",
#                                        split = ",")))
## top centre
#legendPos <- as.numeric(unlist(strsplit("0.38,0.96",
#                                        split = ",")))
## top right
#legendPos <- as.numeric(unlist(strsplit("0.75,0.96",
#                                        split = ",")))
## bottom left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.40",
#                                        split = ",")))
#ChIPNames <- unlist(strsplit("AGO5_D7_sRNA_Rep1,Input_D7_sRNA_Rep1",
#                             split = ","))
#ChIPNamesDir <- unlist(strsplit("miRNAseq_seedling_floral_Bradamante_Gutzat_2022_bioRxiv/snakemake_miRNAseq_,miRNAseq_seedling_floral_Bradamante_Gutzat_2022_bioRxiv/snakemake_miRNAseq_",
#                                split = ","))
#ChIPNamesPlot <- unlist(strsplit("AGO5 D7 Rep1,Input D7 Rep1",
#                                 split = ","))
#ChIPColours <- unlist(strsplit("cyan3,grey50,cyan4,grey40,dodgerblue4,black",
#                               split = ","))
#sRNAsize <- unlist(strsplit("21,22,24",
#                            split = ","))
#ATHILAFam <- "ATHILA1"
#refbase <- "t2t-col.20210610"

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
align <- args[2]
bodyLength <- as.numeric(args[3])
ATHILA_bodyLength <- as.numeric(args[4])
TEsf_bodyLength <- as.numeric(args[5])
upstream <- as.numeric(args[6])
downstream <- as.numeric(args[6])
flankName <- args[7]
binSize <- as.numeric(args[8])
ATHILA_binSize <- as.numeric(args[9])
binName <- args[10]
ATHILA_binName <- args[11]
legendPos <- as.numeric(unlist(strsplit(args[12],
                                        split = ",")))
ChIPNames <- unlist(strsplit(args[13],
                             split = ","))
ChIPNamesDir <- unlist(strsplit(args[14],
                                split = ","))
ChIPNamesPlot <- unlist(strsplit(args[15],
                                 split = ","))
ChIPColours <- unlist(strsplit(args[16],
                               split = ","))
sRNAsize <- unlist(strsplit(args[17],
                            split = ","))
ATHILAFam <- args[18]
refbase <- args[19]

ChIPNames <- rep(ChIPNames, 3)
ChIPNamesDir <- rep(ChIPNamesDir, 3)
sRNAsize <- as.character(sapply(sRNAsize, function(x) rep(x, 2)))
ChIPNamesPlot <- paste0(rep(ChIPNamesPlot, 3), " ", sRNAsize, "-nt")
yLabPlot <- "sRNAs (TPM)"
log2ChIPNames <- ChIPNames
log2ChIPNamesPlot <- ChIPNamesPlot
log2ChIPColours <- ChIPColours

options(stringsAsFactors = F)
library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)
library(extrafontdb)
#extrafont::loadfonts()

outDir <- paste0(paste0(chrName, collapse = "_"), "/")
plotDir <- paste0(outDir, "plots/", ATHILAFam, "/")
plotDir_log2 <- paste0(outDir, "plots/", ATHILAFam, "/log2/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))
system(paste0("[ -d ", plotDir_log2, " ] || mkdir -p ", plotDir_log2))

# Get row indices corresponding to ATHILAFam elements in CENATHILA_BED and nonCENATHILA_BED,
# to be used for extracting rows from coverage matrices
CENATHILA_BED <- read.table(paste0("/home/ajt200/rds/hpc-work/pancentromere/annotation/ATHILA/ATHILA/",
                                   refbase, "/",
                                   "CENATHILA_in_", refbase, "_",
                                   paste0(chrName, collapse = "_"), ".bed"),
                            header = F)
colnames(CENATHILA_BED) <- c("chr", "start0based", "end", "name", "class", "strand")
print(paste0("CENATHILA families in ", refbase, ":"))
#[1] "CENATHILA families in t2t-col.20210610:"
print(sort(unique(CENATHILA_BED$class)))
#[1] "ATHILA1"  "ATHILA2"  "ATHILA5"  "ATHILA6A" "ATHILA6B"
if(ATHILAFam %in% c("ATHILA4")) {
  ATHILAFam_CENATHILA_rowIndices <- which(CENATHILA_BED$class == ATHILAFam)
} else {
  ATHILAFam_CENATHILA_rowIndices <- grep(ATHILAFam, CENATHILA_BED$class)
}

nonCENATHILA_BED <- read.table(paste0("/home/ajt200/rds/hpc-work/pancentromere/annotation/ATHILA/ATHILA/",
                                      refbase, "/",
                                      "nonCENATHILA_in_", refbase, "_",
                                      paste0(chrName, collapse = "_"), ".bed"),
                               header = F)
colnames(nonCENATHILA_BED) <- c("chr", "start0based", "end", "name", "class", "strand")
print(paste0("nonCENATHILA families in ", refbase, ":"))
#[1] "nonCENATHILA families in t2t-col.20210610:"
print(sort(unique(nonCENATHILA_BED$class)))
#[1] "ATHILA0"  "ATHILA1"  "ATHILA2"  "ATHILA3"  "ATHILA4"  "ATHILA4C"
#[7] "ATHILA5"  "ATHILA6A" "ATHILA6B" "ATHILA7A"
if(ATHILAFam %in% c("ATHILA4")) {
  ATHILAFam_nonCENATHILA_rowIndices <- which(nonCENATHILA_BED$class == ATHILAFam)
} else {
  ATHILAFam_nonCENATHILA_rowIndices <- grep(ATHILAFam, nonCENATHILA_BED$class)
}


if(length(ATHILAFam_CENATHILA_rowIndices) > 0) {
  if(length(chrName) == 5) {
    CENATHILANamePlot <- paste0("All CEN", ATHILAFam)
    CENATHILA_CENranLocNamePlot <- paste0("All CENranLoc")
  } else {
    CENATHILANamePlot <- paste0(paste0(chrName, collapse = ","), " CEN", ATHILAFam)
    CENATHILA_CENranLocNamePlot <- paste0(paste0(chrName, collapse = ","), " CENranLoc")
  }
} else {
  if(length(chrName) == 5) {
    CENATHILANamePlot <- paste0("All CENATHILA")
    CENATHILA_CENranLocNamePlot <- paste0("All CENranLoc")
  } else {
    CENATHILANamePlot <- paste0(paste0(chrName, collapse = ","), " CENATHILA")
    CENATHILA_CENranLocNamePlot <- paste0(paste0(chrName, collapse = ","), " CENranLoc")
  }
}

if(length(ATHILAFam_nonCENATHILA_rowIndices) > 0) {
  if(length(chrName) == 5) {
    nonCENATHILANamePlot <- paste0("All nonCEN", ATHILAFam)
    nonCENATHILA_nonCENranLocNamePlot <- paste0("All nonCENranLoc")
  } else {
    nonCENATHILANamePlot <- paste0(paste0(chrName, collapse = ","), " nonCEN", ATHILAFam)
    nonCENATHILA_nonCENranLocNamePlot <- paste0(paste0(chrName, collapse = ","), " nonCENranLoc")
  }
} else {
  if(length(chrName) == 5) {
    nonCENATHILANamePlot <- paste0("All nonCENATHILA")
    nonCENATHILA_nonCENranLocNamePlot <- paste0("All nonCENranLoc")
  } else {
    nonCENATHILANamePlot <- paste0(paste0(chrName, collapse = ","), " nonCENATHILA")
    nonCENATHILA_nonCENranLocNamePlot <- paste0(paste0(chrName, collapse = ","), " nonCENranLoc")
  }
}


# Define feature start and end labels for plotting
featureStartLab <- "Start"
featureEndLab <- "End"
geneStartLab <- "TSS"
geneEndLab <- "TTS"

# Load feature matrices for each chromatin dataset, calculate log2(ChIP/control),
ChIPDirs <- sapply(seq_along(ChIPNames), function(x) {
  paste0("/home/ajt200/rds/hpc-work/",
         ChIPNamesDir[x], refbase,
         "/mapped/")
})

## ChIP
# CENATHILA
ChIP_CENATHILAMats <- mclapply(seq_along(ChIPNames), function(x) {
#  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "ATHILAprofiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_", refbase, "_", align, "_", sRNAsize[x], "nt_sort_norm_CENATHILA_in_",
                                paste0(chrName, collapse = "_"), "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
#  })
}, mc.cores = length(ChIPNames))
# Extract coverage matrix rows that correspond to ATHILAFam elements
if(length(ATHILAFam_CENATHILA_rowIndices) > 0) {
  ChIP_CENATHILAMats <- mclapply(seq_along(ChIP_CENATHILAMats), function(x) {
    ChIP_CENATHILAMats[[x]][c(ATHILAFam_CENATHILA_rowIndices), , drop = FALSE]
  }, mc.cores = length(ChIP_CENATHILAMats))
}

# CENATHILA_CENranLoc
ChIP_CENATHILA_CENranLocMats <- mclapply(seq_along(ChIPNames), function(x) {
#  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "ATHILAprofiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_", refbase, "_", align, "_", sRNAsize[x], "nt_sort_norm_CENATHILA_in_",
                                paste0(chrName, collapse = "_"), "_CENranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
#  })
}, mc.cores = length(ChIPNames))
# Extract coverage matrix rows that correspond to ATHILAFam elements
if(length(ATHILAFam_CENATHILA_rowIndices) > 0) {
  ChIP_CENATHILA_CENranLocMats <- mclapply(seq_along(ChIP_CENATHILA_CENranLocMats), function(x) {
    ChIP_CENATHILA_CENranLocMats[[x]][c(ATHILAFam_CENATHILA_rowIndices), , drop = FALSE]
  }, mc.cores = length(ChIP_CENATHILA_CENranLocMats))
}

# nonCENATHILA
ChIP_nonCENATHILAMats <- mclapply(seq_along(ChIPNames), function(x) {
#  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "ATHILAprofiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_", refbase, "_", align, "_", sRNAsize[x], "nt_sort_norm_nonCENATHILA_in_",
                                paste0(chrName, collapse = "_"), "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
#  })
}, mc.cores = length(ChIPNames))
# Extract coverage matrix rows that correspond to ATHILAFam elements
if(length(ATHILAFam_nonCENATHILA_rowIndices) > 0) {
  ChIP_nonCENATHILAMats <- mclapply(seq_along(ChIP_nonCENATHILAMats), function(x) {
    ChIP_nonCENATHILAMats[[x]][c(ATHILAFam_nonCENATHILA_rowIndices), , drop = FALSE]
  }, mc.cores = length(ChIP_nonCENATHILAMats))
}

# nonCENATHILA_nonCENranLoc
ChIP_nonCENATHILA_nonCENranLocMats <- mclapply(seq_along(ChIPNames), function(x) {
#  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "ATHILAprofiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_", refbase, "_", align, "_", sRNAsize[x], "nt_sort_norm_nonCENATHILA_in_",
                                paste0(chrName, collapse = "_"), "_nonCENranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
#  })
}, mc.cores = length(ChIPNames))
# Extract coverage matrix rows that correspond to ATHILAFam elements
if(length(ATHILAFam_nonCENATHILA_rowIndices) > 0) {
  ChIP_nonCENATHILA_nonCENranLocMats <- mclapply(seq_along(ChIP_nonCENATHILA_nonCENranLocMats), function(x) {
    ChIP_nonCENATHILA_nonCENranLocMats[[x]][c(ATHILAFam_nonCENATHILA_rowIndices), , drop = FALSE]
  }, mc.cores = length(ChIP_nonCENATHILA_nonCENranLocMats))
}


# ChIP
# Add column names
for(x in seq_along(ChIP_CENATHILAMats)) {
  colnames(ChIP_CENATHILAMats[[x]]) <- c(paste0("u", 1:(upstream/ATHILA_binSize)),
                                         paste0("t", ((upstream/ATHILA_binSize)+1):((upstream+ATHILA_bodyLength)/ATHILA_binSize)),
                                         paste0("d", (((upstream+ATHILA_bodyLength)/ATHILA_binSize)+1):(((upstream+ATHILA_bodyLength)/ATHILA_binSize)+(downstream/ATHILA_binSize))))

  colnames(ChIP_CENATHILA_CENranLocMats[[x]]) <- c(paste0("u", 1:(upstream/ATHILA_binSize)),
                                                   paste0("t", ((upstream/ATHILA_binSize)+1):((upstream+ATHILA_bodyLength)/ATHILA_binSize)),
                                                   paste0("d", (((upstream+ATHILA_bodyLength)/ATHILA_binSize)+1):(((upstream+ATHILA_bodyLength)/ATHILA_binSize)+(downstream/ATHILA_binSize))))

  colnames(ChIP_nonCENATHILAMats[[x]]) <- c(paste0("u", 1:(upstream/ATHILA_binSize)),
                                            paste0("t", ((upstream/ATHILA_binSize)+1):((upstream+ATHILA_bodyLength)/ATHILA_binSize)),
                                            paste0("d", (((upstream+ATHILA_bodyLength)/ATHILA_binSize)+1):(((upstream+ATHILA_bodyLength)/ATHILA_binSize)+(downstream/ATHILA_binSize))))

  colnames(ChIP_nonCENATHILA_nonCENranLocMats[[x]]) <- c(paste0("u", 1:(upstream/ATHILA_binSize)),
                                                         paste0("t", ((upstream/ATHILA_binSize)+1):((upstream+ATHILA_bodyLength)/ATHILA_binSize)),
                                                         paste0("d", (((upstream+ATHILA_bodyLength)/ATHILA_binSize)+1):(((upstream+ATHILA_bodyLength)/ATHILA_binSize)+(downstream/ATHILA_binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci 
ChIP_mats <- mclapply(seq_along(ChIP_CENATHILAMats), function(x) {
  list(
       # CENATHILA
       ChIP_CENATHILAMats[[x]],
       # CENATHILA_CENranLoc
       ChIP_CENATHILA_CENranLocMats[[x]],
       # nonCENATHILA
       ChIP_nonCENATHILAMats[[x]],
       # nonCENATHILA_nonCENranLoc
       ChIP_nonCENATHILA_nonCENranLocMats[[x]]
      ) 
}, mc.cores = length(ChIP_CENATHILAMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_ChIP <- mclapply(seq_along(ChIP_mats), function(x) {
  lapply(seq_along(ChIP_mats[[x]]), function(y) {
    data.frame(window = colnames(ChIP_mats[[x]][[y]]),
               t(ChIP_mats[[x]][[y]]))
  })
}, mc.cores = length(ChIP_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_ChIP  <- mclapply(seq_along(wideDFfeature_list_ChIP), function(x) {
  lapply(seq_along(ChIP_mats[[x]]), function(y) {
    gather(data  = wideDFfeature_list_ChIP[[x]][[y]],
           key   = feature,
           value = coverage,
           -window)
  }) 
}, mc.cores = length(wideDFfeature_list_ChIP))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_ChIP)) {
  for(y in seq_along(ChIP_mats[[x]])) {
    tidyDFfeature_list_ChIP[[x]][[y]]$window <- factor(tidyDFfeature_list_ChIP[[x]][[y]]$window,
                                                       levels = as.character(wideDFfeature_list_ChIP[[x]][[y]]$window))
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_ChIP  <- mclapply(seq_along(tidyDFfeature_list_ChIP), function(x) {
  lapply(seq_along(ChIP_mats[[x]]), function(y) {
    data.frame(window = as.character(wideDFfeature_list_ChIP[[x]][[y]]$window),
               n      = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_ChIP[[x]][[y]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_ChIP[[x]][[y]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFfeature_list_ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_ChIP[[x]][[y]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_ChIP))

for(x in seq_along(summaryDFfeature_list_ChIP)) {
  for(y in seq_along(ChIP_mats[[x]])) {
    summaryDFfeature_list_ChIP[[x]][[y]]$window <- factor(summaryDFfeature_list_ChIP[[x]][[y]]$window,
                                                          levels = as.character(wideDFfeature_list_ChIP[[x]][[y]]$window))
    summaryDFfeature_list_ChIP[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_ChIP[[x]][[y]])[1])
    summaryDFfeature_list_ChIP[[x]][[y]]$sem <- summaryDFfeature_list_ChIP[[x]][[y]]$sd/sqrt(summaryDFfeature_list_ChIP[[x]][[y]]$n-1)
    summaryDFfeature_list_ChIP[[x]][[y]]$CI_lower <- summaryDFfeature_list_ChIP[[x]][[y]]$mean -
      qt(0.975, df = summaryDFfeature_list_ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_ChIP[[x]][[y]]$sem
    summaryDFfeature_list_ChIP[[x]][[y]]$CI_upper <- summaryDFfeature_list_ChIP[[x]][[y]]$mean +
      qt(0.975, df = summaryDFfeature_list_ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_ChIP[[x]][[y]]$sem
  }
}

# Convert list of lists summaryDFfeature_list_ChIP into
# a list of single data.frames containing all meta-profiles for plotting
CENATHILATmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[1]]
})
CENATHILA_CENranLocTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[2]]
})
nonCENATHILATmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[3]]
})
nonCENATHILA_nonCENranLocTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[4]]
})
names(CENATHILATmp) <- ChIPNamesPlot
names(CENATHILA_CENranLocTmp) <- ChIPNamesPlot
names(nonCENATHILATmp) <- ChIPNamesPlot
names(nonCENATHILA_nonCENranLocTmp) <- ChIPNamesPlot
summaryDFfeature_ChIP <- list(
  bind_rows(CENATHILATmp, .id = "libName"),
  bind_rows(CENATHILA_CENranLocTmp, .id = "libName"),
  bind_rows(nonCENATHILATmp, .id = "libName"),
  bind_rows(nonCENATHILA_nonCENranLocTmp, .id = "libName")
)  
for(x in seq_along(summaryDFfeature_ChIP)) {
  summaryDFfeature_ChIP[[x]]$libName <- factor(summaryDFfeature_ChIP[[x]]$libName,
                                               levels = ChIPNamesPlot)
}

# Define y-axis limits
#ymin_ChIP <- min(c(summaryDFfeature_ChIP[[1]]$CI_lower,
#                   summaryDFfeature_ChIP[[2]]$CI_lower,
#                   summaryDFfeature_ChIP[[3]]$CI_lower,
#                   summaryDFfeature_ChIP[[4]]$CI_lower),
#                 na.rm = T)
#ymax_ChIP <- max(c(summaryDFfeature_ChIP[[1]]$CI_upper,
#                   summaryDFfeature_ChIP[[2]]$CI_upper,
#                   summaryDFfeature_ChIP[[3]]$CI_upper,
#                   summaryDFfeature_ChIP[[4]]$CI_upper),
#                 na.rm = T)
ymin_ChIP <- min(c(summaryDFfeature_ChIP[[1]]$mean,
                   summaryDFfeature_ChIP[[2]]$mean,
                   summaryDFfeature_ChIP[[3]]$mean,
                   summaryDFfeature_ChIP[[4]]$mean),
                 na.rm = T)
ymax_ChIP <- max(c(summaryDFfeature_ChIP[[1]]$mean,
                   summaryDFfeature_ChIP[[2]]$mean,
                   summaryDFfeature_ChIP[[3]]$mean,
                   summaryDFfeature_ChIP[[4]]$mean),
                 na.rm = T)

# Define legend labels
legendLabs <- lapply(seq_along(ChIPNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(ChIPNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = ChIPColours[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon

## CENATHILA
summaryDFfeature <- summaryDFfeature_ChIP[[1]]
ggObj1_combined_ChIP <- ggplot(data = summaryDFfeature,
                               mapping = aes(x = winNo,
                                             y = mean,
                                             group = libName)
                              ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
#geom_ribbon(data = summaryDFfeature,
#            mapping = aes(ymin = CI_lower,
#                          ymax = CI_upper,
#                          fill = libName),
#            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/ATHILA_binSize)+1,
                            (dim(summaryDFfeature_ChIP[[1]])[1]/length(ChIPNames))-(downstream/ATHILA_binSize),
                            dim(summaryDFfeature_ChIP[[1]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/ATHILA_binSize)+1,
                          (dim(summaryDFfeature_ChIP[[1]])[1]/length(ChIPNames))-(downstream/ATHILA_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote(.(yLabPlot))) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(CENATHILANamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## CENATHILA_CENranLoc
summaryDFfeature <- summaryDFfeature_ChIP[[2]]
ggObj2_combined_ChIP <- ggplot(data = summaryDFfeature,
                               mapping = aes(x = winNo,
                                             y = mean,
                                             group = libName)
                              ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
#geom_ribbon(data = summaryDFfeature,
#            mapping = aes(ymin = CI_lower,
#                          ymax = CI_upper,
#                          fill = libName),
#            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/ATHILA_binSize)+1,
                            (dim(summaryDFfeature_ChIP[[2]])[1]/length(ChIPNames))-(downstream/ATHILA_binSize),
                            dim(summaryDFfeature_ChIP[[2]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/ATHILA_binSize)+1,
                          (dim(summaryDFfeature_ChIP[[2]])[1]/length(ChIPNames))-(downstream/ATHILA_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote(.(yLabPlot))) +
annotation_custom(legendLabs[[1]]) +
annotation_custom(legendLabs[[2]]) +
annotation_custom(legendLabs[[3]]) +
annotation_custom(legendLabs[[4]]) +
annotation_custom(legendLabs[[5]]) +
annotation_custom(legendLabs[[6]]) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(CENATHILA_CENranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## nonCENATHILA
summaryDFfeature <- summaryDFfeature_ChIP[[3]]
ggObj3_combined_ChIP <- ggplot(data = summaryDFfeature,
                               mapping = aes(x = winNo,
                                             y = mean,
                                             group = libName)
                              ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
#geom_ribbon(data = summaryDFfeature,
#            mapping = aes(ymin = CI_lower,
#                          ymax = CI_upper,
#                          fill = libName),
#            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/ATHILA_binSize)+1,
                            (dim(summaryDFfeature_ChIP[[3]])[1]/length(ChIPNames))-(downstream/ATHILA_binSize),
                            dim(summaryDFfeature_ChIP[[3]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/ATHILA_binSize)+1,
                          (dim(summaryDFfeature_ChIP[[3]])[1]/length(ChIPNames))-(downstream/ATHILA_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote(.(yLabPlot))) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(nonCENATHILANamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## nonCENATHILA_nonCENranLoc
summaryDFfeature <- summaryDFfeature_ChIP[[4]]
ggObj4_combined_ChIP <- ggplot(data = summaryDFfeature,
                               mapping = aes(x = winNo,
                                             y = mean,
                                             group = libName)
                              ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
#geom_ribbon(data = summaryDFfeature,
#            mapping = aes(ymin = CI_lower,
#                          ymax = CI_upper,
#                          fill = libName),
#            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/ATHILA_binSize)+1,
                            (dim(summaryDFfeature_ChIP[[4]])[1]/length(ChIPNames))-(downstream/ATHILA_binSize),
                            dim(summaryDFfeature_ChIP[[4]])[1]/length(ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/ATHILA_binSize)+1,
                          (dim(summaryDFfeature_ChIP[[4]])[1]/length(ChIPNames))-(downstream/ATHILA_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote(.(yLabPlot))) +
annotation_custom(legendLabs[[1]]) +
annotation_custom(legendLabs[[2]]) +
annotation_custom(legendLabs[[3]]) +
annotation_custom(legendLabs[[4]]) +
annotation_custom(legendLabs[[5]]) +
annotation_custom(legendLabs[[6]]) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(nonCENATHILA_nonCENranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

ggObjGA_combined <- grid.arrange(grobs = list(
                                              ggObj1_combined_ChIP,
                                              ggObj2_combined_ChIP,
                                              ggObj3_combined_ChIP,
                                              ggObj4_combined_ChIP
                                             ),
                                 layout_matrix = cbind(
                                                       1,
                                                       2,
                                                       3,
                                                       4
                                                      ))
ggsave(paste0(plotDir,
              "sRNA_", paste0(unique(sRNAsize), collapse = "_"), "nt_",
              paste0(unique(ChIPNames), collapse = "_"),
              "_avgProfiles_around",
              "_CEN", ATHILAFam, "Fam",
              "_CENranLoc",
              "_nonCEN", ATHILAFam, "Fam",
              "_nonCENranLoc",
              "_in_", refbase, "_",
              paste0(chrName, collapse = "_"), "_", align, ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 7*length(ggObjGA_combined), limitsize = FALSE)

# log2ChIP
log2ChIP_CENATHILAMats <- mclapply(seq_along(ChIP_CENATHILAMats), function(x) {
  log2((ChIP_CENATHILAMats[[x]]+1))
}, mc.cores = length(ChIP_CENATHILAMats))
log2ChIP_CENATHILA_CENranLocMats <- mclapply(seq_along(ChIP_CENATHILA_CENranLocMats), function(x) {
  log2((ChIP_CENATHILA_CENranLocMats[[x]]+1))
}, mc.cores = length(ChIP_CENATHILA_CENranLocMats))
log2ChIP_nonCENATHILAMats <- mclapply(seq_along(ChIP_nonCENATHILAMats), function(x) {
  log2((ChIP_nonCENATHILAMats[[x]]+1))
}, mc.cores = length(ChIP_nonCENATHILAMats))
log2ChIP_nonCENATHILA_nonCENranLocMats <- mclapply(seq_along(ChIP_nonCENATHILA_nonCENranLocMats), function(x) {
  log2((ChIP_nonCENATHILA_nonCENranLocMats[[x]]+1))
}, mc.cores = length(ChIP_nonCENATHILA_nonCENranLocMats))


# log2ChIP
# Add column names
for(x in seq_along(log2ChIP_CENATHILAMats)) {
  colnames(log2ChIP_CENATHILAMats[[x]]) <- c(paste0("u", 1:(upstream/ATHILA_binSize)),
                                             paste0("t", ((upstream/ATHILA_binSize)+1):((upstream+ATHILA_bodyLength)/ATHILA_binSize)),
                                             paste0("d", (((upstream+ATHILA_bodyLength)/ATHILA_binSize)+1):(((upstream+ATHILA_bodyLength)/ATHILA_binSize)+(downstream/ATHILA_binSize))))

  colnames(log2ChIP_CENATHILA_CENranLocMats[[x]]) <- c(paste0("u", 1:(upstream/ATHILA_binSize)),
                                                       paste0("t", ((upstream/ATHILA_binSize)+1):((upstream+ATHILA_bodyLength)/ATHILA_binSize)),
                                                       paste0("d", (((upstream+ATHILA_bodyLength)/ATHILA_binSize)+1):(((upstream+ATHILA_bodyLength)/ATHILA_binSize)+(downstream/ATHILA_binSize))))

  colnames(log2ChIP_nonCENATHILAMats[[x]]) <- c(paste0("u", 1:(upstream/ATHILA_binSize)),
                                                paste0("t", ((upstream/ATHILA_binSize)+1):((upstream+ATHILA_bodyLength)/ATHILA_binSize)),
                                                paste0("d", (((upstream+ATHILA_bodyLength)/ATHILA_binSize)+1):(((upstream+ATHILA_bodyLength)/ATHILA_binSize)+(downstream/ATHILA_binSize))))

  colnames(log2ChIP_nonCENATHILA_nonCENranLocMats[[x]]) <- c(paste0("u", 1:(upstream/ATHILA_binSize)),
                                                             paste0("t", ((upstream/ATHILA_binSize)+1):((upstream+ATHILA_bodyLength)/ATHILA_binSize)),
                                                             paste0("d", (((upstream+ATHILA_bodyLength)/ATHILA_binSize)+1):(((upstream+ATHILA_bodyLength)/ATHILA_binSize)+(downstream/ATHILA_binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci 
log2ChIP_mats <- mclapply(seq_along(log2ChIP_CENATHILAMats), function(x) {
  list(
       # CENATHILA
       log2ChIP_CENATHILAMats[[x]],
       # CENATHILA_CENranLoc
       log2ChIP_CENATHILA_CENranLocMats[[x]],
       # nonCENATHILA
       log2ChIP_nonCENATHILAMats[[x]],
       # nonCENATHILA_nonCENranLoc
       log2ChIP_nonCENATHILA_nonCENranLocMats[[x]]
      ) 
}, mc.cores = length(log2ChIP_CENATHILAMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_log2ChIP <- mclapply(seq_along(log2ChIP_mats), function(x) {
  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
    data.frame(window = colnames(log2ChIP_mats[[x]][[y]]),
               t(log2ChIP_mats[[x]][[y]]))
  })
}, mc.cores = length(log2ChIP_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_log2ChIP  <- mclapply(seq_along(wideDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
    gather(data  = wideDFfeature_list_log2ChIP[[x]][[y]],
           key   = feature,
           value = coverage,
           -window)
  }) 
}, mc.cores = length(wideDFfeature_list_log2ChIP))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_log2ChIP)) {
  for(y in seq_along(log2ChIP_mats[[x]])) {
    tidyDFfeature_list_log2ChIP[[x]][[y]]$window <- factor(tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                                                       levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window))
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_log2ChIP  <- mclapply(seq_along(tidyDFfeature_list_log2ChIP), function(x) {
  lapply(seq_along(log2ChIP_mats[[x]]), function(y) {
    data.frame(window = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window),
               n      = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFfeature_list_log2ChIP[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_log2ChIP[[x]][[y]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_log2ChIP))

for(x in seq_along(summaryDFfeature_list_log2ChIP)) {
  for(y in seq_along(log2ChIP_mats[[x]])) {
    summaryDFfeature_list_log2ChIP[[x]][[y]]$window <- factor(summaryDFfeature_list_log2ChIP[[x]][[y]]$window,
                                                          levels = as.character(wideDFfeature_list_log2ChIP[[x]][[y]]$window))
    summaryDFfeature_list_log2ChIP[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_log2ChIP[[x]][[y]])[1])
    summaryDFfeature_list_log2ChIP[[x]][[y]]$sem <- summaryDFfeature_list_log2ChIP[[x]][[y]]$sd/sqrt(summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)
    summaryDFfeature_list_log2ChIP[[x]][[y]]$CI_lower <- summaryDFfeature_list_log2ChIP[[x]][[y]]$mean -
      qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]]$sem
    summaryDFfeature_list_log2ChIP[[x]][[y]]$CI_upper <- summaryDFfeature_list_log2ChIP[[x]][[y]]$mean +
      qt(0.975, df = summaryDFfeature_list_log2ChIP[[x]][[y]]$n-1)*summaryDFfeature_list_log2ChIP[[x]][[y]]$sem
  }
}

# Convert list of lists summaryDFfeature_list_log2ChIP into
# a list of single data.frames containing all meta-profiles for plotting
CENATHILATmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[1]]
})
CENATHILA_CENranLocTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[2]]
})
nonCENATHILATmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[3]]
})
nonCENATHILA_nonCENranLocTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[4]]
})
names(CENATHILATmp) <- log2ChIPNamesPlot
names(CENATHILA_CENranLocTmp) <- log2ChIPNamesPlot
names(nonCENATHILATmp) <- log2ChIPNamesPlot
names(nonCENATHILA_nonCENranLocTmp) <- log2ChIPNamesPlot
summaryDFfeature_log2ChIP <- list(
  bind_rows(CENATHILATmp, .id = "libName"),
  bind_rows(CENATHILA_CENranLocTmp, .id = "libName"),
  bind_rows(nonCENATHILATmp, .id = "libName"),
  bind_rows(nonCENATHILA_nonCENranLocTmp, .id = "libName")
)  
for(x in seq_along(summaryDFfeature_log2ChIP)) {
  summaryDFfeature_log2ChIP[[x]]$libName <- factor(summaryDFfeature_log2ChIP[[x]]$libName,
                                               levels = log2ChIPNamesPlot)
}

# Define y-axis limits
#ymin_log2ChIP <- min(c(summaryDFfeature_log2ChIP[[1]]$CI_lower,
#                       summaryDFfeature_log2ChIP[[2]]$CI_lower,
#                       summaryDFfeature_log2ChIP[[3]]$CI_lower,
#                       summaryDFfeature_log2ChIP[[4]]$CI_lower),
#                 na.rm = T)
#ymax_log2ChIP <- max(c(summaryDFfeature_log2ChIP[[1]]$CI_upper,
#                       summaryDFfeature_log2ChIP[[2]]$CI_upper,
#                       summaryDFfeature_log2ChIP[[3]]$CI_upper,
#                       summaryDFfeature_log2ChIP[[4]]$CI_upper),
#                 na.rm = T)
ymin_log2ChIP <- min(c(summaryDFfeature_log2ChIP[[1]]$mean,
                       summaryDFfeature_log2ChIP[[2]]$mean,
                       summaryDFfeature_log2ChIP[[3]]$mean,
                       summaryDFfeature_log2ChIP[[4]]$mean),
                 na.rm = T)
ymax_log2ChIP <- max(c(summaryDFfeature_log2ChIP[[1]]$mean,
                       summaryDFfeature_log2ChIP[[2]]$mean,
                       summaryDFfeature_log2ChIP[[3]]$mean,
                       summaryDFfeature_log2ChIP[[4]]$mean),
                 na.rm = T)

# Define legend labels
legendLabs <- lapply(seq_along(log2ChIPNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(log2ChIPNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = log2ChIPColours[x], fontsize = 18)))
})

# Plot average profiles with 95% CI ribbon

## CENATHILA
summaryDFfeature <- summaryDFfeature_log2ChIP[[1]]
ggObj1_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                               mapping = aes(x = winNo,
                                             y = mean,
                                             group = libName)
                              ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
#geom_ribbon(data = summaryDFfeature,
#            mapping = aes(ymin = CI_lower,
#                          ymax = CI_upper,
#                          fill = libName),
#            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/ATHILA_binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPNames))-(downstream/ATHILA_binSize),
                            dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/ATHILA_binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPNames))-(downstream/ATHILA_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(" * .(yLabPlot) * ")")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(CENATHILANamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## CENATHILA_CENranLoc
summaryDFfeature <- summaryDFfeature_log2ChIP[[2]]
ggObj2_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                               mapping = aes(x = winNo,
                                             y = mean,
                                             group = libName)
                              ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
#geom_ribbon(data = summaryDFfeature,
#            mapping = aes(ymin = CI_lower,
#                          ymax = CI_upper,
#                          fill = libName),
#            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/ATHILA_binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPNames))-(downstream/ATHILA_binSize),
                            dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/ATHILA_binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPNames))-(downstream/ATHILA_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(" * .(yLabPlot) * ")")) +
annotation_custom(legendLabs[[1]]) +
annotation_custom(legendLabs[[2]]) +
annotation_custom(legendLabs[[3]]) +
annotation_custom(legendLabs[[4]]) +
annotation_custom(legendLabs[[5]]) +
annotation_custom(legendLabs[[6]]) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(CENATHILA_CENranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## nonCENATHILA
summaryDFfeature <- summaryDFfeature_log2ChIP[[3]]
ggObj3_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                               mapping = aes(x = winNo,
                                             y = mean,
                                             group = libName)
                              ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
#geom_ribbon(data = summaryDFfeature,
#            mapping = aes(ymin = CI_lower,
#                          ymax = CI_upper,
#                          fill = libName),
#            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/ATHILA_binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[3]])[1]/length(log2ChIPNames))-(downstream/ATHILA_binSize),
                            dim(summaryDFfeature_log2ChIP[[3]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/ATHILA_binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[3]])[1]/length(log2ChIPNames))-(downstream/ATHILA_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(" * .(yLabPlot) * ")")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(nonCENATHILANamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## nonCENATHILA_nonCENranLoc
summaryDFfeature <- summaryDFfeature_log2ChIP[[4]]
ggObj4_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                               mapping = aes(x = winNo,
                                             y = mean,
                                             group = libName)
                              ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
#geom_ribbon(data = summaryDFfeature,
#            mapping = aes(ymin = CI_lower,
#                          ymax = CI_upper,
#                          fill = libName),
#            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%6.3f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/ATHILA_binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[4]])[1]/length(log2ChIPNames))-(downstream/ATHILA_binSize),
                            dim(summaryDFfeature_log2ChIP[[4]])[1]/length(log2ChIPNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/ATHILA_binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[4]])[1]/length(log2ChIPNames))-(downstream/ATHILA_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(" * .(yLabPlot) * ")")) +
annotation_custom(legendLabs[[1]]) +
annotation_custom(legendLabs[[2]]) +
annotation_custom(legendLabs[[3]]) +
annotation_custom(legendLabs[[4]]) +
annotation_custom(legendLabs[[5]]) +
annotation_custom(legendLabs[[6]]) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 3.5, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(nonCENATHILA_nonCENranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

ggObjGA_combined <- grid.arrange(grobs = list(
                                              ggObj1_combined_log2ChIP,
                                              ggObj2_combined_log2ChIP,
                                              ggObj3_combined_log2ChIP,
                                              ggObj4_combined_log2ChIP
                                             ),
                                 layout_matrix = cbind(
                                                       1,
                                                       2,
                                                       3,
                                                       4
                                                      ))
ggsave(paste0(plotDir_log2,
              "sRNA_", paste0(unique(sRNAsize), collapse = "_"), "nt_log2_",
              paste0(unique(log2ChIPNames), collapse = "_"),
              "_avgProfiles_around",
              "_CEN", ATHILAFam, "Fam",
              "_CENranLoc",
              "_nonCEN", ATHILAFam, "Fam",
              "_nonCENranLoc",
              "_in_", refbase, "_",
              paste0(chrName, collapse = "_"), "_", align, ".pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 7*length(ggObjGA_combined), limitsize = FALSE)
