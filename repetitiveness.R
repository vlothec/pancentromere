#!/usr/bin/env Rscript

#to compile C code on the HPC cluster:

# module load gcc/9
# gcc --version
# gcc -o HOR.V3.3 main.c

# module load mafft/7.309

# module load R/4.0.3



# to be run on the repeats output from TRASH. The wrapper will reformat the repeat files and align using mafft so the C script will work properly. Later, the outputs are read back and HOR summary is produced:
# - repeat weighted consensus
# - HOR repetitiveness within species
# - 

#only repeats that have a class assigned are being used

.libPaths("/home/pw457/Rpackages")

library(stringr)
library(base)
library(msa)#
library(Biostrings)
library(seqinr)#
library(doParallel)#




###########
#Functions#
###########

revCompString = function(DNAstr) {return(toupper(toString(reverseComplement(DNAString(DNAstr))))) }

temp.folder =  "/rds/user/pw457/hpc-work/repetitiveness"

cen180.files.directory = "/rds/user/pw457/hpc-work/consensus" 

cen180.files.type = "cen180.consensus.scored" 

HOR.files.directory = "/rds/user/pw457/hpc-work/HOR/single2/horCSVs"

HOR.files.type = "HORs_"

task_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))    #1:330
print(paste("slurm array ID for this job is ", task_id))

setwd(temp.folder)
  
repeats.files = list.files(path = cen180.files.directory, pattern = cen180.files.type, recursive = FALSE) 
  
repeats = read.csv(paste(cen180.files.directory, repeats.files[task_id], sep = "/"), stringsAsFactors = FALSE)

HORfiles = list.files(path = HOR.files.directory, pattern = HOR.files.type, recursive = FALSE) 

found = FALSE
for(i in 1 : length(HORfiles))
{
  HORfiles[i]
  if(grepl(repeats$assembly[1], HORfiles[i]))
  {
    if(grepl(repeats$chromosome[1], HORfiles[i]))
    {
      found = TRUE
      HOR.file.no = i
      break
    }
  }
}

genome = repeats$assembly[1]
chromosome = repeats$chromosome[1]

# names(HORs)
#[1] "X"                        
#[2]"start_A"
#[3] "end_A"                    
#[4]"start_B"
#[5] "end_B"                    
#[6]"direction.1.para_2.perp."
#[7] "total_variant"            
#[8]"start.A.bp"
#[9] "start.B.bp"               
#[10]"end.A.bp"
#[11] "end.B.bp"                 
#[12]"chrA"
#[13] "chrB"                     
#[14]"block.size.in.units"
#[15] "block.A.size.bp"          
#[16]"block.B.size.bp"

#> names(repeats)         #OUTDATED
#[1] "X"                               
#[2]"start"
#[3] "end"                             
#[4]"width"
#[5] "seq"                             
#[6]"strand"
#[7] "class"                           
#[8]"chromosome"
#[9] "est.chromosome.no"               
#[10]"assembly"
#[11] "weighted.consensus.score"        
#[12]"edit.distance"
#[13] "weighted.consensus.outlier"      
#[14]"edit.distance.outlier"
#[15] "weighted.consensus.score.smooth" 
#[16]"edit.distance.smooth"
#[17] "in.region"
if(!file.exists(paste(temp.folder, "/repeats/cen180.consensus.repetitiveness", genome, ".", chromosome, ".csv", sep = "")))
{
  print("running")
  if(found)
  {
    HORs = read.csv(paste(HOR.files.directory, HORfiles[HOR.file.no], sep = "/"), stringsAsFactors = FALSE)
    
    repeats$HORcount = 0
    repeats$HORlengthsSum = 0
    
    for(i in 1 : nrow(HORs))
    {
      for(j in 0 : (HORs$block.size.in.units[i] - 1))
      {
        repeats$HORcount[HORs$start_A[i] + j] = repeats$HORcount[HORs$start_A[i] + j] + 1
        repeats$HORcount[HORs$start_B[i] + j] = repeats$HORcount[HORs$start_B[i] + j] + 1
        repeats$HORlengthsSum[HORs$start_A[i] + j] = repeats$HORlengthsSum[HORs$start_A[i] + j] + HORs$block.size.in.units[i]
        repeats$HORlengthsSum[HORs$start_B[i] + j] = repeats$HORlengthsSum[HORs$start_B[i] + j] + HORs$block.size.in.units[i]
      }
    }
    
    
    
    if(FALSE) #plotting, turn off as regions and outliers not defined
    {
      
      #get regions 
      
      {#regions start
        regions = data.frame(genome = vector(mode = "character", length = 0), chromosome = vector(mode = "character", length = 0), start = vector(mode = "numeric", length = 0), end = vector(mode = "integer", length = 0), region.cen180.count = vector(mode = "integer", length = 0))
        g = 1
        index = 1
        while(g < nrow(repeats))
        {
          region.start = repeats$start[g]
          region.chromosome = repeats$chromosome[g]
          region.genome= repeats$assembly[g]
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
      
      repeats.in.regions = repeats
      genome = repeats$assembly[1]
      chromosome = repeats$chromosome[1]
      
      repeats.in.regions$averagedSUMofHORlengths = 0
      
      slide.window = 25
      
      for(k in 1 : nrow(regions))
      {
        repeats.in.regions.a = repeats.in.regions[repeats.in.regions$in.region == k,]
        repeats.in.regions.b = repeats.in.regions[repeats.in.regions$in.region != k,]
        if(nrow(repeats.in.regions.a) > slide.window)
        {
          for(j in 1 : nrow(repeats.in.regions.a))
          {
            if(j < (slide.window + 1))
            {
              repeats.in.regions.a$averagedSUMofHORlengths[j] = mean(repeats.in.regions.a$HORlengthsSum[1 : (j+slide.window)])
            } else if(j > (nrow(repeats.in.regions) - slide.window - 1))
            {
              repeats.in.regions.a$averagedSUMofHORlengths[j] = mean(repeats.in.regions.a$HORlengthsSum[(j-slide.window) : nrow(repeats.in.regions.a)])
            } else
            {
              repeats.in.regions.a$averagedSUMofHORlengths[j] = mean(repeats.in.regions.a$HORlengthsSum[(j-slide.window) : (j+slide.window)])
            }
          }
        } else
        {
          for(j in 1 : nrow(repeats.in.regions.a))
          {
            repeats.in.regions.a$averagedSUMofHORlengths[j] = mean(repeats.in.regions.a$HORlengthsSum)
          }
        }
        
        repeats.in.regions = rbind(repeats.in.regions.b, repeats.in.regions.a)
        repeats.in.regions = repeats.in.regions[order(repeats.in.regions$start),]
      }
      
      repeats.in.regions$averagedSUMofHORlengths[is.na(repeats.in.regions$averagedSUMofHORlengths)] = repeats.in.regions$HORlengthsSum[is.na(repeats.in.regions$averagedSUMofHORlengths)]
      
      
      png(filename = paste(temp.folder, "/plots/", genome, ".", chromosome, ".png", sep = ""), width = (3000 + nrow(regions)*2000), height = 4500, pointsize = 50)
      par(mfrow = c(3,nrow(regions)*2))
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
               main = paste(genome, chromosome, j, sep = " "))
        } else 
        {
          plot(repeats.in.regions$start[repeats.in.regions$in.region == j], 
               repeats.in.regions$weighted.consensus.score.smooth[repeats.in.regions$in.region == j],
               xlab = "coordinates",
               ylab = "wighted score",
               main = paste(genome, chromosome, j, sep = " "), col = "red")
        } 
        if(nrow(region.repeats.no.outliers) > 100)
        {
          hist(region.repeats.no.outliers$weighted.consensus.score.smooth,
               breaks = "Scott",
               xlab = "wighted score",
               ylab = "frequency",
               main = paste(genome, chromosome, j, sep = " "))
        } else if(length(repeats.in.regions$weighted.consensus.score.smooth[repeats.in.regions$in.region == j]) > 1)
        {
          hist(repeats.in.regions$weighted.consensus.score.smooth[repeats.in.regions$in.region == j],
               breaks = "Scott",
               xlab = "wighted score",
               ylab = "frequency",
               main = paste(genome, chromosome, j, sep = " "))
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
               main = paste(genome, chromosome, j, sep = " "))
        } else 
        {
          plot(repeats.in.regions$start[repeats.in.regions$in.region == j], 
               repeats.in.regions$edit.distance.smooth[repeats.in.regions$in.region == j],
               xlab = "coordinates",
               ylab = "edit distance",
               main = paste(genome, chromosome, j, sep = " "), col = "red")
        }
        if(nrow(region.repeats.no.outliers) > 100)
        {
          hist(region.repeats.no.outliers$edit.distance.smooth,
               breaks = "Scott",
               xlab = "edit distance",
               ylab = "frequency",
               main = paste(genome, chromosome, j, sep = " "))
        } else
        {
          if(length(repeats.in.regions$weighted.consensus.score.smooth[repeats.in.regions$in.region == j]) > 1)
          {
            hist(repeats.in.regions$edit.distance.smooth[repeats.in.regions$in.region == j],
                 breaks = "Scott",
                 xlab = "edit distance",
                 ylab = "frequency",
                 main = paste(genome, chromosome, j, sep = " "))
          } else
          {
            plot.new()
          }
        }
      }
      for(j in 1 : nrow(regions))
      {
        region.repeats.no.outliers = repeats.in.regions[repeats.in.regions$in.region == j,]
        region.repeats.no.outliers = region.repeats.no.outliers[region.repeats.no.outliers$weighted.consensus.outlier == 0,]
        if(nrow(region.repeats.no.outliers) > 100)
        {
          plot(region.repeats.no.outliers$start, 
               region.repeats.no.outliers$averagedSUMofHORlengths,
               type = "l",
               xlab = "coordinates",
               ylab = "sum of HOR lengths",
               main = paste(genome, chromosome, j, sep = " "))
        } else 
        {
          plot(repeats.in.regions$start[repeats.in.regions$in.region == j], 
               repeats.in.regions$averagedSUMofHORlengths[repeats.in.regions$in.region == j],
               xlab = "coordinates",
               ylab = "sum of HOR lengths",
               main = paste(genome, chromosome, j, sep = " "), col = "red")
        } 
        if(nrow(region.repeats.no.outliers) > 100)
        {
          hist(region.repeats.no.outliers$averagedSUMofHORlengths,
               breaks = "Scott",
               xlab = "sum of HOR lengths",
               ylab = "frequency",
               main = paste(genome, chromosome, j, sep = " "))
        } else if(length(repeats.in.regions$averagedSUMofHORlengths[repeats.in.regions$in.region == j]) > 1)
        {
          hist(repeats.in.regions$averagedSUMofHORlengths[repeats.in.regions$in.region == j],
               breaks = "Scott",
               xlab = "sum of HOR lengths",
               ylab = "frequency",
               main = paste(genome, chromosome, j, sep = " "))
        } 
      }
      dev.off()
    }
    
    write.csv(x = repeats, file = paste(temp.folder, "/repeats/cen180.consensus.repetitiveness", genome, ".", chromosome, ".csv", sep = ""))
    #write repeats.in.regions if assigning to regions and plotting
  }
  
} else 
{
  print("skipped")
}


   
print("done")



