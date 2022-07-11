#!/usr/bin/env Rscript


.libPaths("/home/pw457/Rpackages")


library(stringr)
library(base)
library(msa)#
library(Biostrings)
library(seqinr)#
library(stringdist)


revCompString = function(DNAstr) {return(toupper(toString(reverseComplement(DNAString(DNAstr))))) }

temp.folder =  "/rds/user/pw457/hpc-work/consensus"

cen180.files.directory = "/rds/user/pw457/hpc-work/panC_TRASH/repeats_reformat"

setwd(temp.folder)

repeats.files = list.files(path = cen180.files.directory, pattern = "all.repeats.from", recursive = FALSE) 


i = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste("slurm array ID for this job is ", i))



repeats = read.csv(paste(cen180.files.directory, repeats.files[i], sep = "/"), header = TRUE)

repeats$class[is.na(repeats$class)] = 0

repeats = repeats[repeats$class == "aTha178",]

genome = repeats$assembly[1]

chromosomes = unique(repeats$chromosome)

repeatsAll = repeats



for(im in 1 : 1) #for 5 chromosomes separately
{
  if(!file.exists(paste("cen180.consensus.scored", genome, ".", chromosomes[im], ".csv", sep = "")))
  {
    print("doing")
    repeats = repeatsAll[repeatsAll$chromosome == chromosomes[im],]
    
    
    write.fasta(file.out = paste(genome, ".", chromosomes[im], ".repeats.to.align.fasta", sep = ""), 
                names = paste(repeats$chromosome, repeats$start, repeats$end, sep = "_"), sequences = str_split(repeats$sequence_strand_adjusted, pattern = ""), open = "w", as.string = FALSE)
    
    system(paste("/rds/user/pw457/hpc-work/mafft/mafft-linux64/mafft.bat --quiet --retree 2 --inputorder ",
                 paste(genome, ".", chromosomes[im], ".repeats.to.align.fasta", sep = ""), " > ",
                 paste(genome, ".", chromosomes[im], ".aligned.fasta", sep = ""), sep = ""), intern = TRUE)
    
    alignment = read.alignment(paste(genome, ".", chromosomes[im], ".aligned.fasta", sep = ""), format = "fasta")
    
    alignmentVector = vector(mode = "character", length = length(alignment$seq))
    
    for(j in 1: length(alignmentVector))
    {
      alignmentVector[j] = alignment$seq[[j]][1] %>% strsplit(, split = "")
    }
    
    #similarity to consensus based on weighted consensus 
    
    consensusA = vector(mode = "numeric", length = length(alignmentVector[[1]]))
    consensusT = vector(mode = "numeric", length = length(alignmentVector[[1]]))
    consensusC = vector(mode = "numeric", length = length(alignmentVector[[1]]))
    consensusG = vector(mode = "numeric", length = length(alignmentVector[[1]]))
    consensusN = vector(mode = "numeric", length = length(alignmentVector[[1]]))
    
    increment = 1/length(alignmentVector)
    
    for(j in 1 : length(alignmentVector))
    {
      for(k in 1 : length(alignmentVector[[1]]))
      {
        consensusA[k] = consensusA[k] + increment * (alignmentVector[[j]][k] == "a")
        consensusT[k] = consensusT[k] + increment * (alignmentVector[[j]][k] == "t")
        consensusC[k] = consensusC[k] + increment * (alignmentVector[[j]][k] == "c")
        consensusG[k] = consensusG[k] + increment * (alignmentVector[[j]][k] == "g")
        consensusN[k] = consensusN[k] + increment * (alignmentVector[[j]][k] == "-")
      }
    }
    
    cen180.frequencies = data.frame(consensusA, consensusC, consensusG, consensusT)
    cen180.frequencies$sum = cen180.frequencies$consensusA + cen180.frequencies$consensusC + cen180.frequencies$consensusG + cen180.frequencies$consensusT
    cen180.frequencies$consensusN = 1 - cen180.frequencies$sum
    
    repeats$weighted.consensus.score = 0
    for(j in 1 : nrow(repeats))
    {
      for(k in 1 : nrow(cen180.frequencies))
      {
        repeats$weighted.consensus.score[j] = repeats$weighted.consensus.score[j] + cen180.frequencies$consensusA[k]*(alignmentVector[[j]][k] != "a")
        repeats$weighted.consensus.score[j] = repeats$weighted.consensus.score[j] + cen180.frequencies$consensusC[k]*(alignmentVector[[j]][k] != "c")
        repeats$weighted.consensus.score[j] = repeats$weighted.consensus.score[j] + cen180.frequencies$consensusT[k]*(alignmentVector[[j]][k] != "t")
        repeats$weighted.consensus.score[j] = repeats$weighted.consensus.score[j] + cen180.frequencies$consensusG[k]*(alignmentVector[[j]][k] != "g")
        repeats$weighted.consensus.score[j] = repeats$weighted.consensus.score[j] + cen180.frequencies$consensusN[k]*(alignmentVector[[j]][k] != "-")
      }
    }
    
    #similarity to consensus based on edit distance
    
    repeats$edit.distance = 0
    
    consensus180 = consensus(alignment, method = "majority")
    consensus180 = consensus180[consensus180 != "-"]
    consensus180 =toupper(paste(consensus180, collapse = ""))
    
    for(j in 1 : nrow(repeats))
    {
      repeats$edit.distance[j] = stringdist(consensus180, repeats$sequence_strand_adjusted[j], method = "lv")
    }
    
    #mark outliers
    if(FALSE)
    {
      repeats$weighted.consensus.outlier = 0
      repeats$edit.distance.outlier = 0
      
      quantiles.weighted.consensus <- quantile(repeats$weighted.consensus.score, probs=c(.00, .95), na.rm = FALSE)
      quantiles.edit.distance <- quantile(repeats$edit.distance, probs=c(.00, .95), na.rm = FALSE)
      
      repeats$weighted.consensus.outlier[repeats$weighted.consensus.score > quantiles.weighted.consensus[2]] = 1
      repeats$weighted.consensus.outlier[repeats$weighted.consensus.score < quantiles.weighted.consensus[1]] = 1
      repeats$edit.distance.outlier[repeats$edit.distance > quantiles.edit.distance[2]] = 1
      repeats$edit.distance.outlier[repeats$edit.distance > quantiles.edit.distance[2]] = 1
      
      #get regions 
      
      {#regions start
        regions = data.frame(genome = vector(mode = "character", length = 0), chromosome = vector(mode = "character", length = 0), start = vector(mode = "numeric", length = 0), end = vector(mode = "integer", length = 0), region.cen180.count = vector(mode = "integer", length = 0))
        g = 1
        index = 1
        while(g < nrow(repeats))
        {
          region.start = repeats$start[g]
          region.chromosome = repeats$chromosome[g]
          region.genome= repeats$fasta.file.name[g]
          region.end = repeats$end[g]
          region.cen180.count = 0
          g = g + 1
          while((g < nrow(repeats)) && ((repeats$start[g] - repeats$end[g - 1]) < 25000) && (repeats$chromosome[g] == repeats$chromosome[g - 1]))
          {
            region.end = repeats$end[g]
            region.cen180.count = region.cen180.count + 1
            g = g + 1
          }
          region.cen180.count = region.cen180.count + 1
          add = data.frame(region.genome, region.chromosome, region.start, region.end, region.cen180.count)
          names(add) = names(regions)
          regions = rbind(regions, add)
          index = index + 1
        }
        
        regions$size = regions$end - regions$start
        
        regions$distToNext = 0
        for(g in 1 : nrow(regions))
        {
          regions$distToNext[g] = regions$start[g+1] - regions$end[g] 
        }
        
        regions = regions[regions$size > 5000,]
        regions = regions[regions$region.cen180.count > 28,]
        
        if(nrow(regions) > 1)
        {
          for(g in 1 : (nrow(regions) - 1))
          {
            if(nrow(regions) > g)
            {
              if(regions$chromosome[g] != regions$chromosome[g+1])
              {
                regions$distToNext[g] = NA
              }
            }
          }
        }
      }#regions end
      
      
      
      #smooth both
      
      repeats$weighted.consensus.score.smooth = 0
      repeats$edit.distance.smooth = 0
      
      slide.window = 25
      
      repeats.in.regions = repeats
      
      repeats.in.regions$in.region = 0
      
      for(j in nrow(repeats.in.regions) : 1)
      {
        print(j)
        for(k in 1 : nrow(regions))
        {
          if(repeats.in.regions$start[j] > regions$start[k] & repeats.in.regions$start[j] < regions$end[k])
          {
            repeats.in.regions$in.region[j] = k
          }
        }
      }
      
      for(k in 1 : nrow(regions))
      {
        repeats.in.regions.a = repeats.in.regions[repeats.in.regions$in.region == k,]
        repeats.in.regions.b = repeats.in.regions[repeats.in.regions$in.region != k,]
        if(nrow(repeats.in.regions.a) > slide.window)
        {
          for(j in 1 : nrow(repeats.in.regions.a))
          {
            print(j)
            if(j < (slide.window + 1))
            {
              repeats.in.regions.a$weighted.consensus.score.smooth[j] = mean(repeats.in.regions.a$weighted.consensus.score[1 : (j+slide.window)])
              repeats.in.regions.a$edit.distance.smooth[j] = mean(repeats.in.regions.a$edit.distance[1 : (j+slide.window)])
            } else if(j > (nrow(repeats.in.regions) - slide.window - 1))
            {
              repeats.in.regions.a$weighted.consensus.score.smooth[j] = mean(repeats.in.regions.a$weighted.consensus.score[(j-slide.window) : nrow(repeats.in.regions.a)])
              repeats.in.regions.a$edit.distance.smooth[j] = mean(repeats.in.regions.a$edit.distance[(j-slide.window) : nrow(repeats.in.regions.a)])
            } else
            {
              repeats.in.regions.a$weighted.consensus.score.smooth[j] = mean(repeats.in.regions.a$weighted.consensus.score[(j-slide.window) : (j+slide.window)])
              repeats.in.regions.a$edit.distance.smooth[j] = mean(repeats.in.regions.a$edit.distance[(j-slide.window) : (j+slide.window)])
            }
          }
        } else
        {
          for(j in 1 : nrow(repeats.in.regions.a))
          {
            repeats.in.regions.a$weighted.consensus.score.smooth[j] = mean(repeats.in.regions.a$weighted.consensus.score)
            repeats.in.regions.a$edit.distance.smooth[j] = mean(repeats.in.regions.a$edit.distance)
          }
        }
        
        repeats.in.regions = rbind(repeats.in.regions.b, repeats.in.regions.a)
        repeats.in.regions = repeats.in.regions[order(repeats.in.regions$start),]
      }
      
      
      
      repeats.in.regions$weighted.consensus.score.smooth[is.na(repeats.in.regions$weighted.consensus.score.smooth)] = repeats.in.regions$weighted.consensus.score[is.na(repeats.in.regions$weighted.consensus.score.smooth)]
      repeats.in.regions$edit.distance.smooth[is.na(repeats.in.regions$edit.distance.smooth)] = repeats.in.regions$edit.distance[is.na(repeats.in.regions$edit.distance.smooth)]
      
      png(filename = paste("plots/", genome, ".", chromosomes[im], ".png", sep = ""), width = (3000 + nrow(regions)*2000), height = 3000, pointsize = 50)
      par(mfrow = c(2,nrow(regions)*2))
      for(j in 1 : nrow(regions))
      {
        region.repeats.no.outliers = repeats.in.regions[repeats.in.regions$in.region == j,]
        region.repeats.no.outliers = region.repeats.no.outliers[region.repeats.no.outliers$weighted.consensus.outlier == 0,]
        if(nrow(region.repeats.no.outliers) > 100)
        {
          plot(region.repeats.no.outliers$start, 
               region.repeats.no.outliers$weighted.consensus.score.smooth,
               type = "l",
               xlab = "coordinates",
               ylab = "wighted score",
               main = paste(genome, chromosomes[im], j, sep = " "))
        } else 
        {
          plot(repeats.in.regions$start[repeats.in.regions$in.region == j], 
               repeats.in.regions$weighted.consensus.score.smooth[repeats.in.regions$in.region == j],
               xlab = "coordinates",
               ylab = "wighted score",
               main = paste(genome, chromosomes[im], j, sep = " "), col = "red")
        } 
        if(nrow(region.repeats.no.outliers) > 100)
        {
          hist(region.repeats.no.outliers$weighted.consensus.score.smooth,
               breaks = "Scott",
               xlab = "wighted score",
               ylab = "frequency",
               main = paste(genome, chromosomes[im], j, sep = " "))
        } else if(length(repeats.in.regions$weighted.consensus.score.smooth[repeats.in.regions$in.region == j]) > 1)
        {
          hist(repeats.in.regions$weighted.consensus.score.smooth[repeats.in.regions$in.region == j],
               breaks = "Scott",
               xlab = "wighted score",
               ylab = "frequency",
               main = paste(genome, chromosomes[im], j, sep = " "))
        } 
      }
      for(j in 1 : nrow(regions))
      {
        region.repeats.no.outliers = repeats.in.regions[repeats.in.regions$in.region == j,]
        region.repeats.no.outliers = region.repeats.no.outliers[region.repeats.no.outliers$edit.distance.outlier == 0,]
        if(nrow(region.repeats.no.outliers) > 100)
        {
          plot(region.repeats.no.outliers$start, 
               region.repeats.no.outliers$edit.distance.smooth,
               type = "l",
               xlab = "coordinates",
               ylab = "edit distance",
               main = paste(genome, chromosomes[im], j, sep = " "))
        } else 
        {
          plot(repeats.in.regions$start[repeats.in.regions$in.region == j], 
               repeats.in.regions$edit.distance.smooth[repeats.in.regions$in.region == j],
               xlab = "coordinates",
               ylab = "edit distance",
               main = paste(genome, chromosomes[im], j, sep = " "), col = "red")
        }
        if(nrow(region.repeats.no.outliers) > 100)
        {
          hist(region.repeats.no.outliers$edit.distance.smooth,
               breaks = "Scott",
               xlab = "edit distance",
               ylab = "frequency",
               main = paste(genome, chromosomes[im], j, sep = " "))
        } else
        {
          if(length(repeats.in.regions$weighted.consensus.score.smooth[repeats.in.regions$in.region == j]) > 1)
          {
            hist(repeats.in.regions$edit.distance.smooth[repeats.in.regions$in.region == j],
                 breaks = "Scott",
                 xlab = "edit distance",
                 ylab = "frequency",
                 main = paste(genome, chromosomes[im], j, sep = " "))
          } else
          {
            plot.new()
          }
        }
      }
      dev.off()
    }
    
    
    #plot(repeats.in.regions$weighted.consensus.score.smooth, repeats.in.regions$edit.distance.smooth)
    write.csv(x = repeats, file = paste("cen180.consensus.scored", genome, ".", chromosomes[im], ".csv", sep = ""))
    
    
    #remove(repeats.in.regions, regions, alignment, alignmentVector, cen180.frequencies, add, region.repeats.no.outliers, repeats.in.regions.a, repeats.in.regions.b, repeats)
    gc()
  } else
  {
    print("not doing")
  }
}



