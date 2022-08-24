#!/usr/bin/env Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 08.06.2022

# Calculate and plot metaprofiles of ChIP-seq
# (feature windowed means and 95% confidence intervals, CIs)
# for CEN178, ATHILA, Ty3 retrotransposons, and randomly positioned loci


# Usage:
# conda activate R-4.0.0
# ./CEN180_CENranLoc_CENAthila_nonCENAthila_Gypsy_1metaprofile_combined4accessions.R 'Chr1,Chr2,Chr3,Chr4,Chr5' both 180 2000 2000 2000 '2kb' 10 10 10bp 10bp '0.02,0.96' 'Col-0_CENH3_PE150_Rep1_ChIP,Cvi-0_CENH3_PE150_Rep1_ChIP,Ler-0_CENH3_PE150_Rep1_ChIP,Tanz-1_CENH3_PE150_Rep1_ChIP' 'Col-0_CENH3_PE150_Rep1_input,Cvi-0_CENH3_PE150_Rep1_input,Ler-0_CENH3_PE150_Rep1_input,Tanz-1_CENH3_PE150_Rep1_input' 'CENH3_PE150_mn359_20220420/snakemake_ChIPseq_Col-0.ragtag_scaffolds,CENH3_PE150_mn359_20220420/snakemake_ChIPseq_Cvi-0.ragtag_scaffolds,CENH3_PE150_mn359_20220420/snakemake_ChIPseq_Ler-0_110x.ragtag_scaffolds,CENH3_PE150_mn359_20220420/snakemake_ChIPseq_Tanz-1.patch.scaffold.Chr' 'CENH3_PE150_mn359_20220420/snakemake_ChIPseq_Col-0.ragtag_scaffolds,CENH3_PE150_mn359_20220420/snakemake_ChIPseq_Cvi-0.ragtag_scaffolds,CENH3_PE150_mn359_20220420/snakemake_ChIPseq_Ler-0_110x.ragtag_scaffolds,CENH3_PE150_mn359_20220420/snakemake_ChIPseq_Tanz-1.patch.scaffold.Chr' '4 acc. CENH3' '4 acc. input' 'darkorange1' 'darkorange1' 'Col-0.ragtag_scaffolds,Cvi-0.ragtag_scaffolds,Ler-0_110x.ragtag_scaffolds,Tanz-1.patch.scaffold.Chr'
# conda deactivate

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#align <- "both"
#bodyLength <- 180
#Athila_bodyLength <- 2000
#TEsf_bodyLength <- 2000
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#binSize <- 10
#Athila_binSize <- 10
#binName <- "10bp"
#Athila_binName <- "10bp"
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
#ChIPNames <- unlist(strsplit("Col-0_CENH3_PE150_Rep1_ChIP,Cvi-0_CENH3_PE150_Rep1_ChIP,Ler-0_CENH3_PE150_Rep1_ChIP,Tanz-1_CENH3_PE150_Rep1_ChIP",
#                             split = ","))
#controlNames <- unlist(strsplit("Col-0_CENH3_PE150_Rep1_input,Cvi-0_CENH3_PE150_Rep1_input,Ler-0_CENH3_PE150_Rep1_input,Tanz-1_CENH3_PE150_Rep1_input",
#                                split = ","))
#ChIPNamesDir <- unlist(strsplit("CENH3_PE150_mn359_20220420/snakemake_ChIPseq_Col-0.ragtag_scaffolds,CENH3_PE150_mn359_20220420/snakemake_ChIPseq_Cvi-0.ragtag_scaffolds,CENH3_PE150_mn359_20220420/snakemake_ChIPseq_Ler-0_110x.ragtag_scaffolds,CENH3_PE150_mn359_20220420/snakemake_ChIPseq_Tanz-1.patch.scaffold.Chr",
#                                split = ","))
#controlNamesDir <- unlist(strsplit("CENH3_PE150_mn359_20220420/snakemake_ChIPseq_Col-0.ragtag_scaffolds,CENH3_PE150_mn359_20220420/snakemake_ChIPseq_Cvi-0.ragtag_scaffolds,CENH3_PE150_mn359_20220420/snakemake_ChIPseq_Ler-0_110x.ragtag_scaffolds,CENH3_PE150_mn359_20220420/snakemake_ChIPseq_Tanz-1.patch.scaffold.Chr",
#                                   split = ","))
#ChIPNamesPlot <- unlist(strsplit("4 acc. CENH3",
#                                 split = ","))
#controlNamesPlot <- unlist(strsplit("4 acc. input",
#                                    split = ","))
#ChIPColours <- unlist(strsplit("darkorange1",
#                               split = ","))
#controlColours <- unlist(strsplit("darkorange1",
#                                  split = ","))
#refbase <- unlist(strsplit("Col-0.ragtag_scaffolds,Cvi-0.ragtag_scaffolds,Ler-0_110x.ragtag_scaffolds,Tanz-1.patch.scaffold.Chr",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
align <- args[2]
bodyLength <- as.numeric(args[3])
Athila_bodyLength <- as.numeric(args[4])
TEsf_bodyLength <- as.numeric(args[5])
upstream <- as.numeric(args[6])
downstream <- as.numeric(args[6])
flankName <- args[7]
binSize <- as.numeric(args[8])
Athila_binSize <- as.numeric(args[9])
binName <- args[10]
Athila_binName <- args[11]
legendPos <- as.numeric(unlist(strsplit(args[12],
                                        split = ",")))
ChIPNames <- unlist(strsplit(args[13],
                             split = ","))
controlNames <- unlist(strsplit(args[14],
                                split = ","))
ChIPNamesDir <- unlist(strsplit(args[15],
                                split = ","))
controlNamesDir <- unlist(strsplit(args[16],
                                   split = ","))
ChIPNamesPlot <- unlist(strsplit(args[17],
                                 split = ","))
controlNamesPlot <- unlist(strsplit(args[18],
                                    split = ","))
ChIPColours <- unlist(strsplit(args[19],
                               split = ","))
controlColours <- unlist(strsplit(args[20],
                                  split = ","))
refbase <- unlist(strsplit(args[21],
                           split = ","))

options(stringsAsFactors = F)
library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)
#extrafont::loadfonts()

outDir <- paste0(paste0(chrName, collapse = "_"), "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir -p ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir -p ", plotDir))

if(length(chrName) == 5) {
  featureNamePlot <- "All CEN178"
  ranLocNamePlot <- "All CENranLoc"
  AthilaNamePlot <- "All CENATHILA"
  nonCENAthilaNamePlot <- "All nonCENATHILA"
  TEsfNamePlot <- "All Ty3"
} else {
  featureNamePlot <- paste0(paste0(chrName, collapse = ","), " CEN178")
  ranLocNamePlot <- paste0(paste0(chrName, collapse = ","), " CENranLoc")
  AthilaNamePlot <- paste0(paste0(chrName, collapse = ","), " CENATHILA")
  nonCENAthilaNamePlot <- paste0(paste0(chrName, collapse = ","), " nonCENATHILA")
  TEsfNamePlot <- paste0(paste0(chrName, collapse = ","), " Ty3")
}

# Define feature start and end labels for plotting
featureStartLab <- "Start"
featureEndLab <- "End"
geneStartLab <- "TSS"
geneEndLab <- "TTS"

log2ChIPNames <- ChIPNames
log2ChIPNamesPlot <- ChIPNamesPlot
log2ChIPColours <- ChIPColours

# Load feature matrices for each chromatin dataset, calculate log2(ChIP/control),
ChIPDirs <- sapply(seq_along(ChIPNames), function(x) {
  paste0("/home/ajt200/analysis/",
         ChIPNamesDir[x],
         "/mapped/")
})

controlDirs <- sapply(seq_along(controlNames), function(x) {
  paste0("/home/ajt200/analysis/",
         controlNamesDir[x],
         "/mapped/")
})

## ChIP
# feature
ChIP_featureMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "CEN180profiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_CEN180_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3,
                         colClasses = c(rep("NULL", 1000/binSize), rep(NA, ((1000*2)+(180))/binSize), rep("NULL", 1000/binSize))))
  })
}, mc.cores = length(ChIPNames))
# If features from multiple chromosomes are to be analysed,
# concatenate the corresponding feature coverage matrices
ChIP_featureMats <- mclapply(seq_along(ChIP_featureMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_featureMats[[x]])
  } else {
    ChIP_featureMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_featureMats))

## control
# feature
control_featureMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x], "CEN180profiles/matrices/",
                                controlNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_CEN180_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3,
                         colClasses = c(rep("NULL", 1000/binSize), rep(NA, ((1000*2)+(180))/binSize), rep("NULL", 1000/binSize))))
  })
}, mc.cores = length(controlNames))
# If features from multiple chromosomes are to be analysed,
# concatenate the corresponding feature coverage matrices
control_featureMats <- mclapply(seq_along(control_featureMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_featureMats[[x]])
  } else {
    control_featureMats[[x]][[1]]
  }
}, mc.cores = length(control_featureMats))


# Due to variable presence/absence of ATHILA per chromosome,
# load full matrix across all chromosomes

## ChIP
# ranLoc
ChIP_ranLocMats <- mclapply(seq_along(ChIPNames), function(x) {
  as.matrix(read.table(paste0(ChIPDirs[x], "ATHILAprofiles/matrices/",
                              ChIPNames[x],
                              "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_CENATHILA_in_",
                              paste0(chrName, collapse = "_"), "_CENranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames))

## control
# ranLoc
control_ranLocMats <- mclapply(seq_along(controlNames), function(x) {
  as.matrix(read.table(paste0(controlDirs[x], "ATHILAprofiles/matrices/",
                              controlNames[x],
                              "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_CENATHILA_in_",
                              paste0(chrName, collapse = "_"), "_CENranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(controlNames))


# Due to variable presence/absence of ATHILA per chromosome,
# load full matrix across all chromosomes

## ChIP
# Athila
ChIP_AthilaMats <- mclapply(seq_along(ChIPNames), function(x) {
  as.matrix(read.table(paste0(ChIPDirs[x], "ATHILAprofiles/matrices/",
                              ChIPNames[x],
                              "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_CENATHILA_in_",
                              paste0(chrName, collapse = "_"), "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames))

## control
# Athila
control_AthilaMats <- mclapply(seq_along(controlNames), function(x) {
  as.matrix(read.table(paste0(controlDirs[x], "ATHILAprofiles/matrices/",
                              controlNames[x],
                              "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_CENATHILA_in_",
                              paste0(chrName, collapse = "_"), "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(controlNames))


## ChIP
# nonCENAthila
ChIP_nonCENAthilaMats <- mclapply(seq_along(ChIPNames), function(x) {
  as.matrix(read.table(paste0(ChIPDirs[x], "ATHILAprofiles/matrices/",
                              ChIPNames[x],
                              "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_nonCENATHILA_in_",
                              paste0(chrName, collapse = "_"), "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(ChIPNames))

## control
# nonCENAthila
control_nonCENAthilaMats <- mclapply(seq_along(controlNames), function(x) {
  as.matrix(read.table(paste0(controlDirs[x], "ATHILAprofiles/matrices/",
                              controlNames[x],
                              "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_nonCENATHILA_in_",
                              paste0(chrName, collapse = "_"), "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                       header = F, skip = 3))
}, mc.cores = length(controlNames))


## ChIP
# TEsf
ChIP_TEsfMats <- mclapply(seq_along(ChIPNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(ChIPDirs[x], "TEprofiles/matrices/",
                                ChIPNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_TEs_Gypsy_LTR_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(ChIPNames))
# If TEsfs from multiple chromosomes are to be analysed,
# concatenate the corresponding TEsf coverage matrices
ChIP_TEsfMats <- mclapply(seq_along(ChIP_TEsfMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, ChIP_TEsfMats[[x]])
  } else {
    ChIP_TEsfMats[[x]][[1]]
  }
}, mc.cores = length(ChIP_TEsfMats))

## control
# TEsf
control_TEsfMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    as.matrix(read.table(paste0(controlDirs[x], "TEprofiles/matrices/",
                                controlNames[x],
                                "_MappedOn_", refbase[x], "_lowXM_", align, "_sort_norm_TEs_Gypsy_LTR_in_",
                                chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  })
}, mc.cores = length(controlNames))
# If TEsfs from multiple chromosomes are to be analysed,
# concatenate the corresponding TEsf coverage matrices
control_TEsfMats <- mclapply(seq_along(control_TEsfMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_TEsfMats[[x]])
  } else {
    control_TEsfMats[[x]][[1]]
  }
}, mc.cores = length(control_TEsfMats))

# Conditionally calculate log2(ChIP/control)
# for each matrix depending on library
# feature
log2ChIP_featureMats <- mclapply(seq_along(ChIP_featureMats), function(x) {
  print(paste0(ChIPNames[x], " library; using ", controlNames[x], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_featureMats[[x]]+1)/(control_featureMats[[x]]+1))
}, mc.cores = length(ChIP_featureMats))
log2ChIP_featureMats <- list(do.call(rbind, log2ChIP_featureMats))
ChIP_featureMats <- list(do.call(rbind, ChIP_featureMats))
control_featureMats <- list(do.call(rbind, control_featureMats))

# Conditionally calculate log2(ChIP/control)
# for each matrix depending on library
# ranLoc
log2ChIP_ranLocMats <- mclapply(seq_along(ChIP_ranLocMats), function(x) {
  print(paste0(ChIPNames[x], " library; using ", controlNames[x], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_ranLocMats[[x]]+1)/(control_ranLocMats[[x]]+1))
}, mc.cores = length(ChIP_ranLocMats))
log2ChIP_ranLocMats <- list(do.call(rbind, log2ChIP_ranLocMats))
ChIP_ranLocMats <- list(do.call(rbind, ChIP_ranLocMats))
control_ranLocMats <- list(do.call(rbind, control_ranLocMats))

# Conditionally calculate log2(ChIP/control)
# for each matrix depending on library
# Athila
log2ChIP_AthilaMats <- mclapply(seq_along(ChIP_AthilaMats), function(x) {
  print(paste0(ChIPNames[x], " library; using ", controlNames[x], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_AthilaMats[[x]]+1)/(control_AthilaMats[[x]]+1))
}, mc.cores = length(ChIP_AthilaMats))
log2ChIP_AthilaMats <- list(do.call(rbind, log2ChIP_AthilaMats))
ChIP_AthilaMats <- list(do.call(rbind, ChIP_AthilaMats))
control_AthilaMats <- list(do.call(rbind, control_AthilaMats))

# Conditionally calculate log2(ChIP/control)
# for each matrix depending on library
# nonCENAthila
log2ChIP_nonCENAthilaMats <- mclapply(seq_along(ChIP_nonCENAthilaMats), function(x) {
  print(paste0(ChIPNames[x], " library; using ", controlNames[x], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_nonCENAthilaMats[[x]]+1)/(control_nonCENAthilaMats[[x]]+1))
}, mc.cores = length(ChIP_nonCENAthilaMats))
log2ChIP_nonCENAthilaMats <- list(do.call(rbind, log2ChIP_nonCENAthilaMats))
ChIP_nonCENAthilaMats <- list(do.call(rbind, ChIP_nonCENAthilaMats))
control_nonCENAthilaMats <- list(do.call(rbind, control_nonCENAthilaMats))

# Conditionally calculate log2(ChIP/control)
# for each matrix depending on library
# TEsf
log2ChIP_TEsfMats <- mclapply(seq_along(ChIP_TEsfMats), function(x) {
  print(paste0(ChIPNames[x], " library; using ", controlNames[x], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_TEsfMats[[x]]+1)/(control_TEsfMats[[x]]+1))
}, mc.cores = length(ChIP_TEsfMats))
log2ChIP_TEsfMats <- list(do.call(rbind, log2ChIP_TEsfMats))
ChIP_TEsfMats <- list(do.call(rbind, ChIP_TEsfMats))
control_TEsfMats <- list(do.call(rbind, control_TEsfMats))

# log2ChIP
# Add column names
for(x in seq_along(log2ChIP_featureMats)) {
  colnames(log2ChIP_featureMats[[x]]) <- c(paste0("u", 1:((upstream-1000)/binSize)),
                                           paste0("t", (((upstream-1000)/binSize)+1):(((upstream-1000)+bodyLength)/binSize)),
                                           paste0("d", ((((upstream-1000)+bodyLength)/binSize)+1):((((upstream-1000)+bodyLength)/binSize)+((downstream-1000)/binSize))))
  colnames(log2ChIP_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/Athila_binSize)),
                                          paste0("t", ((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)),
                                          paste0("d", (((upstream+Athila_bodyLength)/Athila_binSize)+1):(((upstream+Athila_bodyLength)/Athila_binSize)+(downstream/Athila_binSize))))
  colnames(log2ChIP_AthilaMats[[x]]) <- c(paste0("u", 1:(upstream/Athila_binSize)),
                                          paste0("t", ((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)),
                                          paste0("d", (((upstream+Athila_bodyLength)/Athila_binSize)+1):(((upstream+Athila_bodyLength)/Athila_binSize)+(downstream/Athila_binSize))))
  colnames(log2ChIP_nonCENAthilaMats[[x]]) <- c(paste0("u", 1:(upstream/Athila_binSize)),
                                                paste0("t", ((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)),
                                                paste0("d", (((upstream+Athila_bodyLength)/Athila_binSize)+1):(((upstream+Athila_bodyLength)/Athila_binSize)+(downstream/Athila_binSize))))
  colnames(log2ChIP_TEsfMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                        paste0("t", ((upstream/binSize)+1):((upstream+TEsf_bodyLength)/binSize)),
                                        paste0("d", (((upstream+TEsf_bodyLength)/binSize)+1):(((upstream+TEsf_bodyLength)/binSize)+(downstream/binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci
log2ChIP_mats <- mclapply(seq_along(log2ChIP_featureMats), function(x) {
  list(
       # features
       log2ChIP_featureMats[[x]],
       # ranLocs
       log2ChIP_ranLocMats[[x]],
       # Athilas
       log2ChIP_AthilaMats[[x]],
       # nonCENAthilas
       log2ChIP_nonCENAthilaMats[[x]],
       # TEsfs
       log2ChIP_TEsfMats[[x]]
      )
}, mc.cores = length(log2ChIP_featureMats))

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
featureTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[1]]
})
ranLocTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[2]]
})
AthilaTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[3]]
})
nonCENAthilaTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[4]]
})
TEsfTmp <- lapply(seq_along(summaryDFfeature_list_log2ChIP), function(x) {
  summaryDFfeature_list_log2ChIP[[x]][[5]]
})
names(featureTmp) <- log2ChIPNamesPlot
names(ranLocTmp) <- log2ChIPNamesPlot
names(AthilaTmp) <- log2ChIPNamesPlot
names(nonCENAthilaTmp) <- log2ChIPNamesPlot
names(TEsfTmp) <- log2ChIPNamesPlot
summaryDFfeature_log2ChIP <- list(
  bind_rows(featureTmp, .id = "libName"),
  bind_rows(ranLocTmp, .id = "libName"),
  bind_rows(AthilaTmp, .id = "libName"),
  bind_rows(nonCENAthilaTmp, .id = "libName"),
  bind_rows(TEsfTmp, .id = "libName")
)
for(x in seq_along(summaryDFfeature_log2ChIP)) {
  summaryDFfeature_log2ChIP[[x]]$libName <- factor(summaryDFfeature_log2ChIP[[x]]$libName,
                                                   levels = log2ChIPNamesPlot)
}

# Define y-axis limits
ymin_log2ChIP <- min(c(summaryDFfeature_log2ChIP[[1]]$CI_lower,
                       summaryDFfeature_log2ChIP[[2]]$CI_lower,
                       summaryDFfeature_log2ChIP[[3]]$CI_lower,
                       summaryDFfeature_log2ChIP[[4]]$CI_lower,
                       summaryDFfeature_log2ChIP[[5]]$CI_lower),
                 na.rm = T)
ymax_log2ChIP <- max(c(summaryDFfeature_log2ChIP[[1]]$CI_upper,
                       summaryDFfeature_log2ChIP[[2]]$CI_upper,
                       summaryDFfeature_log2ChIP[[3]]$CI_upper,
                       summaryDFfeature_log2ChIP[[4]]$CI_upper,
                       summaryDFfeature_log2ChIP[[5]]$CI_upper),
                 na.rm = T)
#ymin_log2ChIP <- -0.1
#ymax_log2ChIP <- 2.8

log2ChIPPlaceholder <- "None"
ChIPPlaceholder <- "None"
controlPlaceholder <- "None"

# Define legend labels
legendLabs <- lapply(seq_along(log2ChIPNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(log2ChIPNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = log2ChIPColours[x], fontsize = 22)))
})

# Plot average profiles with 95% CI ribbon
## feature
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
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%3.1f", x)) +
scale_x_discrete(breaks = c(1,
                            ((upstream-1000)/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPPlaceholder))-((downstream-1000)/binSize),
                            dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPPlaceholder)),
                 labels = c(paste0("-", "1kb"),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", "1kb"))) +
geom_vline(xintercept = c(((upstream-1000)/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[1]])[1]/length(log2ChIPPlaceholder))-((downstream-1000)/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(featureNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## ranLoc
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
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%3.1f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/Athila_binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPPlaceholder))-(downstream/Athila_binSize),
                            dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPPlaceholder)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/Athila_binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[2]])[1]/length(log2ChIPPlaceholder))-(downstream/Athila_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
annotation_custom(legendLabs[[1]]) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Athila
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
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%3.1f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/Athila_binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[3]])[1]/length(log2ChIPPlaceholder))-(downstream/Athila_binSize),
                            dim(summaryDFfeature_log2ChIP[[3]])[1]/length(log2ChIPPlaceholder)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/Athila_binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[3]])[1]/length(log2ChIPPlaceholder))-(downstream/Athila_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(AthilaNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## nonCENAthila
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
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%3.1f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/Athila_binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[4]])[1]/length(log2ChIPPlaceholder))-(downstream/Athila_binSize),
                            dim(summaryDFfeature_log2ChIP[[4]])[1]/length(log2ChIPPlaceholder)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/Athila_binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[4]])[1]/length(log2ChIPPlaceholder))-(downstream/Athila_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(nonCENAthilaNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## TEsf
summaryDFfeature <- summaryDFfeature_log2ChIP[[5]]
ggObj5_combined_log2ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = log2ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = log2ChIPColours) +
scale_y_continuous(limits = c(ymin_log2ChIP, ymax_log2ChIP),
                   labels = function(x) sprintf("%3.1f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_log2ChIP[[5]])[1]/length(log2ChIPPlaceholder))-(downstream/binSize),
                            dim(summaryDFfeature_log2ChIP[[5]])[1]/length(log2ChIPPlaceholder)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_log2ChIP[[5]])[1]/length(log2ChIPPlaceholder))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(TEsfNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

ggObjGA_combined <- grid.arrange(grobs = list(
                                              ggObj1_combined_log2ChIP,
                                              ggObj2_combined_log2ChIP,
                                              ggObj3_combined_log2ChIP,
                                              ggObj4_combined_log2ChIP,
                                              ggObj5_combined_log2ChIP
                                             ),
                                 layout_matrix = cbind(
                                                       1,
                                                       2,
                                                       3,
                                                       4,
                                                       5
                                                      ))
ggsave(paste0(plotDir,
              "log2ChIPcontrol_",
              paste0(log2ChIPNames, collapse = "_"),
              "_avgProfiles_around",
              "_CEN180_CENranLoc_CENAthila_nonCENAthila_Gypsy_",
              paste0(chrName, collapse = "_"), "_", align, "_combined.pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 6.5*5, limitsize = FALSE)


# ChIP
# Add column names
for(x in seq_along(ChIP_featureMats)) {
  colnames(ChIP_featureMats[[x]]) <- c(paste0("u", 1:((upstream-1000)/binSize)),
                                       paste0("t", (((upstream-1000)/binSize)+1):(((upstream-1000)+bodyLength)/binSize)),
                                       paste0("d", ((((upstream-1000)+bodyLength)/binSize)+1):((((upstream-1000)+bodyLength)/binSize)+((downstream-1000)/binSize))))
  colnames(ChIP_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/Athila_binSize)),
                                      paste0("t", ((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)),
                                      paste0("d", (((upstream+Athila_bodyLength)/Athila_binSize)+1):(((upstream+Athila_bodyLength)/Athila_binSize)+(downstream/Athila_binSize))))
  colnames(ChIP_AthilaMats[[x]]) <- c(paste0("u", 1:(upstream/Athila_binSize)),
                                      paste0("t", ((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)),
                                      paste0("d", (((upstream+Athila_bodyLength)/Athila_binSize)+1):(((upstream+Athila_bodyLength)/Athila_binSize)+(downstream/Athila_binSize))))
  colnames(ChIP_nonCENAthilaMats[[x]]) <- c(paste0("u", 1:(upstream/Athila_binSize)),
                                      paste0("t", ((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)),
                                      paste0("d", (((upstream+Athila_bodyLength)/Athila_binSize)+1):(((upstream+Athila_bodyLength)/Athila_binSize)+(downstream/Athila_binSize))))
  colnames(ChIP_TEsfMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                    paste0("t", ((upstream/binSize)+1):((upstream+TEsf_bodyLength)/binSize)),
                                    paste0("d", (((upstream+TEsf_bodyLength)/binSize)+1):(((upstream+TEsf_bodyLength)/binSize)+(downstream/binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci
ChIP_mats <- mclapply(seq_along(ChIP_featureMats), function(x) {
  list(
       # features
       ChIP_featureMats[[x]],
       # ranLocs
       ChIP_ranLocMats[[x]],
       # Athilas
       ChIP_AthilaMats[[x]],
       # nonCENAthilas
       ChIP_nonCENAthilaMats[[x]],
       # TEsfs
       ChIP_TEsfMats[[x]]
      )
}, mc.cores = length(ChIP_featureMats))

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
featureTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[1]]
})
ranLocTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[2]]
})
AthilaTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[3]]
})
nonCENAthilaTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[4]]
})
TEsfTmp <- lapply(seq_along(summaryDFfeature_list_ChIP), function(x) {
  summaryDFfeature_list_ChIP[[x]][[5]]
})
names(featureTmp) <- ChIPNamesPlot
names(ranLocTmp) <- ChIPNamesPlot
names(AthilaTmp) <- ChIPNamesPlot
names(nonCENAthilaTmp) <- ChIPNamesPlot
names(TEsfTmp) <- ChIPNamesPlot
summaryDFfeature_ChIP <- list(
  bind_rows(featureTmp, .id = "libName"),
  bind_rows(ranLocTmp, .id = "libName"),
  bind_rows(AthilaTmp, .id = "libName"),
  bind_rows(nonCENAthilaTmp, .id = "libName"),
  bind_rows(TEsfTmp, .id = "libName")
)
for(x in seq_along(summaryDFfeature_ChIP)) {
  summaryDFfeature_ChIP[[x]]$libName <- factor(summaryDFfeature_ChIP[[x]]$libName,
                                               levels = ChIPNamesPlot)
}

# Define y-axis limits
ymin_ChIP <- min(c(summaryDFfeature_ChIP[[1]]$CI_lower,
                       summaryDFfeature_ChIP[[2]]$CI_lower,
                       summaryDFfeature_ChIP[[3]]$CI_lower,
                       summaryDFfeature_ChIP[[4]]$CI_lower,
                       summaryDFfeature_ChIP[[5]]$CI_lower),
                 na.rm = T)
ymax_ChIP <- max(c(summaryDFfeature_ChIP[[1]]$CI_upper,
                       summaryDFfeature_ChIP[[2]]$CI_upper,
                       summaryDFfeature_ChIP[[3]]$CI_upper,
                       summaryDFfeature_ChIP[[4]]$CI_upper,
                       summaryDFfeature_ChIP[[5]]$CI_upper),
                 na.rm = T)
#ymin_ChIP <- -1.0
#ymax_ChIP <- 29.0

# Define legend labels
legendLabs <- lapply(seq_along(ChIPNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(ChIPNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = ChIPColours[x], fontsize = 22)))
})

# Plot average profiles with 95% CI ribbon
## feature
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
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%3.1f", x)) +
scale_x_discrete(breaks = c(1,
                            ((upstream-1000)/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[1]])[1]/length(ChIPPlaceholder))-((downstream-1000)/binSize),
                            dim(summaryDFfeature_ChIP[[1]])[1]/length(ChIPPlaceholder)),
                 labels = c(paste0("-", "1kb"),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", "1kb"))) +
geom_vline(xintercept = c(((upstream-1000)/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[1]])[1]/length(ChIPPlaceholder))-((downstream-1000)/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Log"[2] * "(ChIP/input)")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(featureNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## ranLoc
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
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%3.1f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/Athila_binSize)+1,
                            (dim(summaryDFfeature_ChIP[[2]])[1]/length(ChIPPlaceholder))-(downstream/Athila_binSize),
                            dim(summaryDFfeature_ChIP[[2]])[1]/length(ChIPPlaceholder)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/Athila_binSize)+1,
                          (dim(summaryDFfeature_ChIP[[2]])[1]/length(ChIPPlaceholder))-(downstream/Athila_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("ChIP")) +
annotation_custom(legendLabs[[1]]) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Athila
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
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%3.1f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/Athila_binSize)+1,
                            (dim(summaryDFfeature_ChIP[[3]])[1]/length(ChIPPlaceholder))-(downstream/Athila_binSize),
                            dim(summaryDFfeature_ChIP[[3]])[1]/length(ChIPPlaceholder)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/Athila_binSize)+1,
                          (dim(summaryDFfeature_ChIP[[3]])[1]/length(ChIPPlaceholder))-(downstream/Athila_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("ChIP")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(AthilaNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## nonCENAthila
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
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%3.1f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/Athila_binSize)+1,
                            (dim(summaryDFfeature_ChIP[[4]])[1]/length(ChIPPlaceholder))-(downstream/Athila_binSize),
                            dim(summaryDFfeature_ChIP[[4]])[1]/length(ChIPPlaceholder)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/Athila_binSize)+1,
                          (dim(summaryDFfeature_ChIP[[4]])[1]/length(ChIPPlaceholder))-(downstream/Athila_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("ChIP")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(nonCENAthilaNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## TEsf
summaryDFfeature <- summaryDFfeature_ChIP[[5]]
ggObj5_combined_ChIP <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = ChIPColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = ChIPColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%3.1f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[5]])[1]/length(ChIPPlaceholder))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[5]])[1]/length(ChIPPlaceholder)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[5]])[1]/length(ChIPPlaceholder))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("ChIP")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(TEsfNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

ggObjGA_combined <- grid.arrange(grobs = list(
                                              ggObj1_combined_ChIP,
                                              ggObj2_combined_ChIP,
                                              ggObj3_combined_ChIP,
                                              ggObj4_combined_ChIP,
                                              ggObj5_combined_ChIP
                                             ),
                                 layout_matrix = cbind(
                                                       1,
                                                       2,
                                                       3,
                                                       4,
                                                       5
                                                      ))
ggsave(paste0(plotDir,
              "ChIP_",
              paste0(ChIPNames, collapse = "_"),
              "_avgProfiles_around",
              "_CEN180_CENranLoc_CENAthila_nonCENAthila_Gypsy_",
              paste0(chrName, collapse = "_"), "_", align, "_combined.pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 6.5*5, limitsize = FALSE)


# control
# Add column names
for(x in seq_along(control_featureMats)) {
  colnames(control_featureMats[[x]]) <- c(paste0("u", 1:((upstream-1000)/binSize)),
                                       paste0("t", (((upstream-1000)/binSize)+1):(((upstream-1000)+bodyLength)/binSize)),
                                       paste0("d", ((((upstream-1000)+bodyLength)/binSize)+1):((((upstream-1000)+bodyLength)/binSize)+((downstream-1000)/binSize))))
  colnames(control_ranLocMats[[x]]) <- c(paste0("u", 1:(upstream/Athila_binSize)),
                                      paste0("t", ((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)),
                                      paste0("d", (((upstream+Athila_bodyLength)/Athila_binSize)+1):(((upstream+Athila_bodyLength)/Athila_binSize)+(downstream/Athila_binSize))))
  colnames(control_AthilaMats[[x]]) <- c(paste0("u", 1:(upstream/Athila_binSize)),
                                      paste0("t", ((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)),
                                      paste0("d", (((upstream+Athila_bodyLength)/Athila_binSize)+1):(((upstream+Athila_bodyLength)/Athila_binSize)+(downstream/Athila_binSize))))
  colnames(control_nonCENAthilaMats[[x]]) <- c(paste0("u", 1:(upstream/Athila_binSize)),
                                      paste0("t", ((upstream/Athila_binSize)+1):((upstream+Athila_bodyLength)/Athila_binSize)),
                                      paste0("d", (((upstream+Athila_bodyLength)/Athila_binSize)+1):(((upstream+Athila_bodyLength)/Athila_binSize)+(downstream/Athila_binSize))))
  colnames(control_TEsfMats[[x]]) <- c(paste0("u", 1:(upstream/binSize)),
                                    paste0("t", ((upstream/binSize)+1):((upstream+TEsf_bodyLength)/binSize)),
                                    paste0("d", (((upstream+TEsf_bodyLength)/binSize)+1):(((upstream+TEsf_bodyLength)/binSize)+(downstream/binSize))))
}

# Create list of lists in which each element in the enclosing list corresponds to a library
# and the two elements in the nested list correspond to coverage matrices for features and random loci
control_mats <- mclapply(seq_along(control_featureMats), function(x) {
  list(
       # features
       control_featureMats[[x]],
       # ranLocs
       control_ranLocMats[[x]],
       # Athilas
       control_AthilaMats[[x]],
       # nonCENAthilas
       control_nonCENAthilaMats[[x]],
       # TEsfs
       control_TEsfMats[[x]]
      )
}, mc.cores = length(control_featureMats))

# Transpose matrix and convert into dataframe
# in which first column is window name
wideDFfeature_list_control <- mclapply(seq_along(control_mats), function(x) {
  lapply(seq_along(control_mats[[x]]), function(y) {
    data.frame(window = colnames(control_mats[[x]][[y]]),
               t(control_mats[[x]][[y]]))
  })
}, mc.cores = length(control_mats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list_control  <- mclapply(seq_along(wideDFfeature_list_control), function(x) {
  lapply(seq_along(control_mats[[x]]), function(y) {
    gather(data  = wideDFfeature_list_control[[x]][[y]],
           key   = feature,
           value = coverage,
           -window)
  }) 
}, mc.cores = length(wideDFfeature_list_control))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list_control)) {
  for(y in seq_along(control_mats[[x]])) {
    tidyDFfeature_list_control[[x]][[y]]$window <- factor(tidyDFfeature_list_control[[x]][[y]]$window,
                                                           levels = as.character(wideDFfeature_list_control[[x]][[y]]$window))
  }
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list_control  <- mclapply(seq_along(tidyDFfeature_list_control), function(x) {
  lapply(seq_along(control_mats[[x]]), function(y) {
    data.frame(window = as.character(wideDFfeature_list_control[[x]][[y]]$window),
               n      = tapply(X     = tidyDFfeature_list_control[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_control[[x]][[y]]$window,
                               FUN   = length),
               mean   = tapply(X     = tidyDFfeature_list_control[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_control[[x]][[y]]$window,
                               FUN   = mean,
                               na.rm = TRUE),
               sd     = tapply(X     = tidyDFfeature_list_control[[x]][[y]]$coverage,
                               INDEX = tidyDFfeature_list_control[[x]][[y]]$window,
                               FUN   = sd,
                               na.rm = TRUE))
  })
}, mc.cores = length(tidyDFfeature_list_control))

for(x in seq_along(summaryDFfeature_list_control)) {
  for(y in seq_along(control_mats[[x]])) {
    summaryDFfeature_list_control[[x]][[y]]$window <- factor(summaryDFfeature_list_control[[x]][[y]]$window,
                                                             levels = as.character(wideDFfeature_list_control[[x]][[y]]$window))
    summaryDFfeature_list_control[[x]][[y]]$winNo <- factor(1:dim(summaryDFfeature_list_control[[x]][[y]])[1])
    summaryDFfeature_list_control[[x]][[y]]$sem <- summaryDFfeature_list_control[[x]][[y]]$sd/sqrt(summaryDFfeature_list_control[[x]][[y]]$n-1)
    summaryDFfeature_list_control[[x]][[y]]$CI_lower <- summaryDFfeature_list_control[[x]][[y]]$mean -
      qt(0.975, df = summaryDFfeature_list_control[[x]][[y]]$n-1)*summaryDFfeature_list_control[[x]][[y]]$sem
    summaryDFfeature_list_control[[x]][[y]]$CI_upper <- summaryDFfeature_list_control[[x]][[y]]$mean +
      qt(0.975, df = summaryDFfeature_list_control[[x]][[y]]$n-1)*summaryDFfeature_list_control[[x]][[y]]$sem
  }
}

# Convert list of lists summaryDFfeature_list_control into
# a list of single data.frames containing all meta-profiles for plotting
featureTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[1]]
})
ranLocTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[2]]
})
AthilaTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[3]]
})
nonCENAthilaTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[4]]
})
TEsfTmp <- lapply(seq_along(summaryDFfeature_list_control), function(x) {
  summaryDFfeature_list_control[[x]][[5]]
})
names(featureTmp) <- controlNamesPlot
names(ranLocTmp) <- controlNamesPlot
names(AthilaTmp) <- controlNamesPlot
names(nonCENAthilaTmp) <- controlNamesPlot
names(TEsfTmp) <- controlNamesPlot
summaryDFfeature_control <- list(
  bind_rows(featureTmp, .id = "libName"),
  bind_rows(ranLocTmp, .id = "libName"),
  bind_rows(AthilaTmp, .id = "libName"),
  bind_rows(nonCENAthilaTmp, .id = "libName"),
  bind_rows(TEsfTmp, .id = "libName")
)
for(x in seq_along(summaryDFfeature_control)) {
  summaryDFfeature_control[[x]]$libName <- factor(summaryDFfeature_control[[x]]$libName,
                                                   levels = controlNamesPlot)
}

# Define y-axis limits
ymin_control <- min(c(summaryDFfeature_control[[1]]$CI_lower,
                      summaryDFfeature_control[[2]]$CI_lower,
                      summaryDFfeature_control[[3]]$CI_lower,
                      summaryDFfeature_control[[4]]$CI_lower,
                      summaryDFfeature_control[[5]]$CI_lower),
                    na.rm = T)
ymax_control <- max(c(summaryDFfeature_control[[1]]$CI_upper,
                      summaryDFfeature_control[[2]]$CI_upper,
                      summaryDFfeature_control[[3]]$CI_upper,
                      summaryDFfeature_control[[4]]$CI_upper,
                      summaryDFfeature_control[[5]]$CI_upper),
                    na.rm = T)
#ymin_control <- -1.0
#ymax_control <- 29.0

# Define legend labels
legendLabs <- lapply(seq_along(controlNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(controlNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = controlColours[x], fontsize = 22)))
})

# Plot average profiles with 95% CI ribbon
## feature
summaryDFfeature <- summaryDFfeature_control[[1]]
ggObj1_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%3.1f", x)) +
scale_x_discrete(breaks = c(1,
                            ((upstream-1000)/binSize)+1,
                            (dim(summaryDFfeature_control[[1]])[1]/length(controlPlaceholder))-((downstream-1000)/binSize),
                            dim(summaryDFfeature_control[[1]])[1]/length(controlPlaceholder)),
                 labels = c(paste0("-", "1kb"),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", "1kb"))) +
geom_vline(xintercept = c(((upstream-1000)/binSize)+1,
                          (dim(summaryDFfeature_control[[1]])[1]/length(controlPlaceholder))-((downstream-1000)/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Control")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(featureNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## ranLoc
summaryDFfeature <- summaryDFfeature_control[[2]]
ggObj2_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%3.1f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/Athila_binSize)+1,
                            (dim(summaryDFfeature_control[[2]])[1]/length(controlPlaceholder))-(downstream/Athila_binSize),
                            dim(summaryDFfeature_control[[2]])[1]/length(controlPlaceholder)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/Athila_binSize)+1,
                          (dim(summaryDFfeature_control[[2]])[1]/length(controlPlaceholder))-(downstream/Athila_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Control")) +
annotation_custom(legendLabs[[1]]) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(ranLocNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## Athila
summaryDFfeature <- summaryDFfeature_control[[3]]
ggObj3_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%3.1f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/Athila_binSize)+1,
                            (dim(summaryDFfeature_control[[3]])[1]/length(controlPlaceholder))-(downstream/Athila_binSize),
                            dim(summaryDFfeature_control[[3]])[1]/length(controlPlaceholder)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/Athila_binSize)+1,
                          (dim(summaryDFfeature_control[[3]])[1]/length(controlPlaceholder))-(downstream/Athila_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Control")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(AthilaNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## nonCENAthila
summaryDFfeature <- summaryDFfeature_control[[4]]
ggObj4_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%3.1f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/Athila_binSize)+1,
                            (dim(summaryDFfeature_control[[4]])[1]/length(controlPlaceholder))-(downstream/Athila_binSize),
                            dim(summaryDFfeature_control[[4]])[1]/length(controlPlaceholder)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/Athila_binSize)+1,
                          (dim(summaryDFfeature_control[[4]])[1]/length(controlPlaceholder))-(downstream/Athila_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Control")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(nonCENAthilaNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

## TEsf
summaryDFfeature <- summaryDFfeature_control[[5]]
ggObj5_combined_control <- ggplot(data = summaryDFfeature,
                                   mapping = aes(x = winNo,
                                                 y = mean,
                                                 group = libName)
                                  ) +
geom_line(data = summaryDFfeature,
          mapping = aes(colour = libName),
          size = 1) +
scale_colour_manual(values = controlColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = controlColours) +
scale_y_continuous(limits = c(ymin_control, ymax_control),
                   labels = function(x) sprintf("%3.1f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_control[[5]])[1]/length(controlPlaceholder))-(downstream/binSize),
                            dim(summaryDFfeature_control[[5]])[1]/length(controlPlaceholder)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_control[[5]])[1]/length(controlPlaceholder))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("Control")) +
theme_bw() +
theme(
      axis.ticks = element_line(size = 1.0, colour = "black"),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.x = element_text(size = 22, colour = "black"),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title = element_text(size = 30, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.0, colour = "black"),
      panel.background = element_blank(),
      plot.margin = unit(c(0.3,1.2,0.0,0.3), "cm"),
      plot.title = element_text(hjust = 0.5, size = 30)) +
ggtitle(bquote(.(TEsfNamePlot) ~ "(" * italic("n") ~ "=" ~
               .(prettyNum(summaryDFfeature$n[1],
                           big.mark = ",", trim = T)) *
               ")"))

ggObjGA_combined <- grid.arrange(grobs = list(
                                              ggObj1_combined_control,
                                              ggObj2_combined_control,
                                              ggObj3_combined_control,
                                              ggObj4_combined_control,
                                              ggObj5_combined_control
                                             ),
                                 layout_matrix = cbind(
                                                       1,
                                                       2,
                                                       3,
                                                       4,
                                                       5
                                                      ))
ggsave(paste0(plotDir,
              "control_",
              paste0(controlNames, collapse = "_"),
              "_avgProfiles_around",
              "_CEN180_CENranLoc_CENAthila_nonCENAthila_Gypsy_",
              paste0(chrName, collapse = "_"), "_", align, "_combined.pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 6.5*5, limitsize = FALSE)
