#!/usr/bin/env Rscript


.libPaths("/home/pw457/Rpackages")

library(stringr)
library(base)
library(msa)#
library(Biostrings)
library(seqinr)#
library(doParallel)#
library(ggplot2)


chromosome = 1

revCompString = function(DNAstr) {return(toupper(toString(reverseComplement(DNAString(DNAstr))))) }



setwd("/rds/user/pw457/hpc-work/HORs_between_chr_across_all")

HOR.files = list.files(path = "/rds/user/pw457/hpc-work/data/HOR_relative_temp", full.names = TRUE)


chromosome_windows = vector(mode = "numeric", length = 100)






for(i in 1 : length(HOR.files))
{
  print(i)
  repeatsAsave = read.csv(file = HOR.files[i])
  if(length(which(is.infinite(repeatsAsave$HORrelativeSharedToTotalA))) > length(which(is.infinite(repeatsAsave$HORrelativeSharedToTotalB))))
  {
    repeatsAsave$HORrelativeSharedToTotalA = repeatsAsave$HORrelativeSharedToTotalB
  }
  
  repeatsAsave$HORrelativeSharedToTotalA[is.infinite(repeatsAsave$HORrelativeSharedToTotalA)] = 0
  
  for(j in 1 : length(chromosome_windows))
  {
    windowA = round(nrow(repeatsAsave)/length(chromosome_windows), digits = 0)
    startA = (j - 1) * windowA + 1
    endA = j * windowA
    if(startA < 1)
    {
      startA = 1
    }
    if(endA > nrow(repeatsAsave))
    {
      endA = nrow(repeatsAsave)
    }
    chromosome_windows[j] = chromosome_windows[j] + mean(repeatsAsave$HORrelativeSharedToTotalA[startA:endA])
  }
  
}

chromosome_windows = chromosome_windows / length(HOR.files)

pdf(file = "plot_test.pdf")
plot(chromosome_windows, type = "l", 
     main = "Syntenic regions across all chr1", 
     xlab = "Chromosome bins (100)", 
     ylab = "HOR based synteny score (0:1)")
dev.off()
































