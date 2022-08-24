#!/usr/bin/env Rscript

# author: Andy Tock
# contact: ajt200@cam.ac.uk
# date: 08.06.2022

# Calculate and plot metaprofiles of DNA methylation proportions
# (CEN180 windowed means and 95% confidence intervals, CIs)
# for CEN178, ATHILA, Ty3 retrotransposons, and randomly positioned loci

# Usage:
# conda activate R-4.0.0
# ./CEN180_CENranLoc_CENAthila_nonCENAthila_Gypsy_DNAmethAllContexts_3metaprofiles_combined4accessions.R 'Chr1,Chr2,Chr3,Chr4,Chr5' 180 2000 2000 2000 2kb 10 10 10bp 10bp '0.02,0.62' 'Col-0_deepsignalDNAmeth_20kb_MappedOn_Col-0.ragtag_scaffolds,Cvi-0_deepsignalDNAmeth_10kb_MappedOn_Cvi-0.ragtag_scaffolds,Ler-0_deepsignalDNAmeth_10kb_MappedOn_Ler-0_110x.ragtag_scaffolds,Tanz-1_deepsignalDNAmeth_unkb_MappedOn_Tanz-1.patch.scaffold.Chr' 'pancentromere/deepsignal_DNAmeth/Col-0.ragtag_scaffolds,pancentromere/deepsignal_DNAmeth/Cvi-0.ragtag_scaffolds,pancentromere/deepsignal_DNAmeth/Ler-0_110x.ragtag_scaffolds,pancentromere/deepsignal_DNAmeth/Tanz-1.patch.scaffold.Chr' 'CpG,CHG,CHH' 'dodgerblue4,dodgerblue1,cyan2' 'Col-0.ragtag_scaffolds,Cvi-0.ragtag_scaffolds,Ler-0_110x.ragtag_scaffolds,Tanz-1.patch.scaffold.Chr'
# conda deactivate

#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
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
## centre left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.62",
#                                        split = ",")))
## bottom left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.40",
#                                        split = ",")))
#ChIPNames <- unlist(strsplit("Col-0_deepsignalDNAmeth_20kb_MappedOn_Col-0.ragtag_scaffolds,Cvi-0_deepsignalDNAmeth_10kb_MappedOn_Cvi-0.ragtag_scaffolds,Ler-0_deepsignalDNAmeth_10kb_MappedOn_Ler-0_110x.ragtag_scaffolds,Tanz-1_deepsignalDNAmeth_unkb_MappedOn_Tanz-1.patch.scaffold.Chr",
#                             split = ","))
#ChIPNamesDir <- unlist(strsplit("pancentromere/deepsignal_DNAmeth/Col-0.ragtag_scaffolds,pancentromere/deepsignal_DNAmeth/Cvi-0.ragtag_scaffolds,pancentromere/deepsignal_DNAmeth/Ler-0_110x.ragtag_scaffolds,pancentromere/deepsignal_DNAmeth/Tanz-1.patch.scaffold.Chr",
#                                split = ","))
#contextNames <- unlist(strsplit("CpG,CHG,CHH",
#                                split = ","))
#contextColours <- unlist(strsplit("dodgerblue4,dodgerblue1,cyan2",
#                               split = ","))
#refbase <- unlist(strsplit("Col-0.ragtag_scaffolds,Cvi-0.ragtag_scaffolds,Ler-0_110x.ragtag_scaffolds,Tanz-1.patch.scaffold.Chr",
#                           split = ","))

args <- commandArgs(trailingOnly = T)
chrName <- unlist(strsplit(args[1],
                           split = ","))
bodyLength <- as.numeric(args[2])
Athila_bodyLength <- as.numeric(args[3])
TEsf_bodyLength <- as.numeric(args[4])
upstream <- as.numeric(args[5])
downstream <- as.numeric(args[5])
flankName <- args[6]
binSize <- as.numeric(args[7])
Athila_binSize <- as.numeric(args[8])
binName <- args[9]
Athila_binName <- args[10]
legendPos <- as.numeric(unlist(strsplit(args[11],
                                        split = ",")))
ChIPNames <- unlist(strsplit(args[12],
                             split = ","))
ChIPNamesDir <- unlist(strsplit(args[13],
                                split = ","))
contextNames <- unlist(strsplit(args[14],
                                split = ","))
contextColours <- unlist(strsplit(args[15],
                                  split = ","))
refbase <- unlist(strsplit(args[16],
                           split = ","))

contextNamesPlot <- gsub("p", "", contextNames)

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

# Load feature matrices for each dataset
ChIPDirs <- sapply(seq_along(ChIPNamesDir), function(x) {
  if(grepl("deepsignal", ChIPNamesDir[x])) {
    paste0("/home/ajt200/analysis/",
           ChIPNamesDir[x],
           "/")
  } else {
    paste0("/home/ajt200/analysis/",
           ChIPNamesDir[x],
           "/coverage/")
  }
})

## ChIP
# feature
ChIP_featureMats <- lapply(seq_along(contextNames), function(w) {
  mclapply(seq_along(ChIPNames), function(x) {
    lapply(seq_along(chrName), function(y) {
      as.matrix(read.table(paste0(ChIPDirs[x], "CEN180profiles/matrices/",
                                  ChIPNames[x], "_", contextNames[w],
                                  "_CEN180_in_",
                                  chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                           header = F, skip = 3,
                           colClasses = c(rep("NULL", 1000/binSize), rep(NA, ((1000*2)+(180))/binSize), rep("NULL", 1000/binSize))))
    })
  }, mc.cores = length(ChIPNames))
})
# If features from multiple chromosomes are to be analysed,
# concatenate the corresponding feature coverage matrices
ChIP_featureMats <- lapply(seq_along(ChIP_featureMats), function(w) {
  tmp <- mclapply(seq_along(ChIP_featureMats[[w]]), function(x) {
    if(length(chrName) > 1) {
      do.call(rbind, ChIP_featureMats[[w]][[x]])
    } else {
      ChIP_featureMats[[w]][[x]][[1]]
    }
  }, mc.cores = length(ChIP_featureMats[[w]]))
  do.call(rbind, tmp)
})

## ChIP
# ranLoc
# Due to variable presence/absence of ATHILA per chromosome,
# load full matrix across all chromosomes
ChIP_ranLocMats <- lapply(seq_along(contextNames), function(w) {
  mclapply(seq_along(ChIPNames), function(x) {
    as.matrix(read.table(paste0(ChIPDirs[x], "ATHILAprofiles/matrices/",
                                ChIPNames[x], "_", contextNames[w],
                                "_CENATHILA_in_",
                                paste0(chrName, collapse = "_"), "_CENranLoc_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  }, mc.cores = length(ChIPNames))
})
ChIP_ranLocMats <- lapply(seq_along(ChIP_ranLocMats), function(w) {
  do.call(rbind, ChIP_ranLocMats[[w]])
})

# ChIP
# Athila
# Due to variable presence/absence of ATHILA per chromosome,
# load full matrix across all chromosomes
ChIP_AthilaMats <- lapply(seq_along(contextNames), function(w) {
  mclapply(seq_along(ChIPNames), function(x) {
    as.matrix(read.table(paste0(ChIPDirs[x], "ATHILAprofiles/matrices/",
                                ChIPNames[x], "_", contextNames[w],
                                "_CENATHILA_in_",
                                paste0(chrName, collapse = "_"), "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  }, mc.cores = length(ChIPNames))
})
ChIP_AthilaMats <- lapply(seq_along(ChIP_AthilaMats), function(w) {
  do.call(rbind, ChIP_AthilaMats[[w]])
})

# ChIP
# nonCENAthila
# Due to variable presence/absence of ATHILA per chromosome,
# load full matrix across all chromosomes
ChIP_nonCENAthilaMats <- lapply(seq_along(contextNames), function(w) {
  mclapply(seq_along(ChIPNames), function(x) {
    as.matrix(read.table(paste0(ChIPDirs[x], "ATHILAprofiles/matrices/",
                                ChIPNames[x], "_", contextNames[w],
                                "_nonCENATHILA_in_",
                                paste0(chrName, collapse = "_"), "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                         header = F, skip = 3))
  }, mc.cores = length(ChIPNames))
})
ChIP_nonCENAthilaMats <- lapply(seq_along(ChIP_nonCENAthilaMats), function(w) {
  do.call(rbind, ChIP_nonCENAthilaMats[[w]])
})

## ChIP
# TEsf
ChIP_TEsfMats <- lapply(seq_along(contextNames), function(w) {
  mclapply(seq_along(ChIPNames), function(x) {
    lapply(seq_along(chrName), function(y) {
      as.matrix(read.table(paste0(ChIPDirs[x], "TEprofiles/matrices/",
                                  ChIPNames[x], "_", contextNames[w],
                                  "_TEs_Gypsy_LTR_in_",
                                  chrName[y], "_matrix_bin", binSize, "bp_flank", flankName, ".tab"),
                           header = F, skip = 3))
    })
  }, mc.cores = length(ChIPNames))
})
# If TEsfs from multiple chromosomes are to be analysed,
# concatenate the corresponding TEsf coverage matrices
ChIP_TEsfMats <- lapply(seq_along(ChIP_TEsfMats), function(w) {
  tmp <- mclapply(seq_along(ChIP_TEsfMats[[w]]), function(x) {
    if(length(chrName) > 1) {
      do.call(rbind, ChIP_TEsfMats[[w]][[x]])
    } else {
      ChIP_TEsfMats[[w]][[x]][[1]]
    }
  }, mc.cores = length(ChIP_TEsfMats[[w]]))
  do.call(rbind, tmp)
})

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
names(featureTmp) <- contextNamesPlot
names(ranLocTmp) <- contextNamesPlot
names(AthilaTmp) <- contextNamesPlot
names(nonCENAthilaTmp) <- contextNamesPlot
names(TEsfTmp) <- contextNamesPlot
summaryDFfeature_ChIP <- list(
  bind_rows(featureTmp, .id = "libName"),
  bind_rows(ranLocTmp, .id = "libName"),
  bind_rows(AthilaTmp, .id = "libName"),
  bind_rows(nonCENAthilaTmp, .id = "libName"),
  bind_rows(TEsfTmp, .id = "libName")
)
for(x in seq_along(summaryDFfeature_ChIP)) {
  summaryDFfeature_ChIP[[x]]$libName <- factor(summaryDFfeature_ChIP[[x]]$libName,
                                               levels = contextNamesPlot)
}

# Define y-axis limits
#ymin_ChIP <- min(c(summaryDFfeature_ChIP[[1]]$CI_lower,
#                   summaryDFfeature_ChIP[[2]]$CI_lower,
#                   summaryDFfeature_ChIP[[3]]$CI_lower,
#                   summaryDFfeature_ChIP[[4]]$CI_lower,
#                   summaryDFfeature_ChIP[[5]]$CI_lower),
#                 na.rm = T)
#ymax_ChIP <- max(c(summaryDFfeature_ChIP[[1]]$CI_upper,
#                   summaryDFfeature_ChIP[[2]]$CI_upper,
#                   summaryDFfeature_ChIP[[3]]$CI_upper,
#                   summaryDFfeature_ChIP[[4]]$CI_upper,
#                   summaryDFfeature_ChIP[[5]]$CI_upper),
#                 na.rm = T)
ymin_ChIP <- 0
ymax_ChIP <- 100

# Define legend labels
legendLabs <- lapply(seq_along(contextNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(contextNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.06), just = "left",
                    gp = gpar(col = contextColours[x], fontsize = 22)))
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
scale_colour_manual(values = contextColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = contextColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%3.0f", x)) +
scale_x_discrete(breaks = c(1,
                            ((upstream-1000)/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[1]])[1]/length(contextNames))-((downstream-1000)/binSize),
                            dim(summaryDFfeature_ChIP[[1]])[1]/length(contextNames)),
                 labels = c(paste0("-", "1kb"),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", "1kb"))) +
geom_vline(xintercept = c(((upstream-1000)/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[1]])[1]/length(contextNames))-((downstream-1000)/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("DNA methylation (%)")) +
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
scale_colour_manual(values = contextColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = contextColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%3.0f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/Athila_binSize)+1,
                            (dim(summaryDFfeature_ChIP[[2]])[1]/length(contextNames))-(downstream/Athila_binSize),
                            dim(summaryDFfeature_ChIP[[2]])[1]/length(contextNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/Athila_binSize)+1,
                          (dim(summaryDFfeature_ChIP[[2]])[1]/length(contextNames))-(downstream/Athila_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("DNA methylation (%)")) +
annotation_custom(legendLabs[[1]]) +
annotation_custom(legendLabs[[2]]) +
annotation_custom(legendLabs[[3]]) +
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
scale_colour_manual(values = contextColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = contextColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%3.0f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/Athila_binSize)+1,
                            (dim(summaryDFfeature_ChIP[[3]])[1]/length(contextNames))-(downstream/Athila_binSize),
                            dim(summaryDFfeature_ChIP[[3]])[1]/length(contextNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/Athila_binSize)+1,
                          (dim(summaryDFfeature_ChIP[[3]])[1]/length(contextNames))-(downstream/Athila_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("DNA methylation (%)")) +
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
scale_colour_manual(values = contextColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = contextColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%3.0f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/Athila_binSize)+1,
                            (dim(summaryDFfeature_ChIP[[4]])[1]/length(contextNames))-(downstream/Athila_binSize),
                            dim(summaryDFfeature_ChIP[[4]])[1]/length(contextNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/Athila_binSize)+1,
                          (dim(summaryDFfeature_ChIP[[4]])[1]/length(contextNames))-(downstream/Athila_binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("DNA methylation (%)")) +
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
scale_colour_manual(values = contextColours) +
geom_ribbon(data = summaryDFfeature,
            mapping = aes(ymin = CI_lower,
                          ymax = CI_upper,
                          fill = libName),
            alpha = 0.4) +
scale_fill_manual(values = contextColours) +
scale_y_continuous(limits = c(ymin_ChIP, ymax_ChIP),
                   labels = function(x) sprintf("%3.0f", x)) +
scale_x_discrete(breaks = c(1,
                            (upstream/binSize)+1,
                            (dim(summaryDFfeature_ChIP[[5]])[1]/length(contextNames))-(downstream/binSize),
                            dim(summaryDFfeature_ChIP[[5]])[1]/length(contextNames)),
                 labels = c(paste0("-", flankName),
                            featureStartLab,
                            featureEndLab,
                            paste0("+", flankName))) +
geom_vline(xintercept = c((upstream/binSize)+1,
                          (dim(summaryDFfeature_ChIP[[5]])[1]/length(contextNames))-(downstream/binSize)),
           linetype = "dashed",
           size = 1) +
labs(x = "",
     y = bquote("DNA methylation (%)")) +
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
              "DNAmethAllContexts_",
              paste0(gsub("_deepsignalDNAmeth_.+", "", ChIPNames), collapse = "_"),
              "_avgProfiles_around",
              "_CEN180_CENranLoc_CENAthila_nonCENAthila_Gypsy_in_",
              paste0(chrName, collapse = "_"), "_combined.pdf"),
       plot = ggObjGA_combined,
       height = 6.5, width = 6.5*5, limitsize = FALSE)
