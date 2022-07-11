#!/usr/bin/env Rscript


.libPaths("/home/pw457/Rpackages")


library(stringr)
library(base)
library(msa)#
library(Biostrings)
library(seqinr)#
library(stringdist)


revCompString = function(DNAstr) {return(toupper(toString(reverseComplement(DNAString(DNAstr))))) }


temp.folder =  "/rds/user/pw457/hpc-work/variant_coordinates"

alignments.files.directory = "/rds/user/pw457/hpc-work/consensusVariants/alignments"

cen180.files.directory = "/rds/user/pw457/hpc-work/panC_TRASH"

setwd(temp.folder)

assemblies.directory = ""

repeats.files = list.files(path = cen180.files.directory, pattern = "all.repeats.from", recursive = FALSE) 

ID = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste("slurm array ID for this job is ", ID))

repeats = read.csv(paste(cen180.files.directory, repeats.files[ID], sep = "/"), header = TRUE)

repeats = repeats[!is.na(repeats$class),]
repeats = repeats[repeats$class == "aTha178",]

genome = strsplit(repeats.files[ID], split = "all.repeats.from.")[[1]][2]
genome = strsplit(genome, split = ".csv")[[1]][1]

output = paste("/rds/user/pw457/hpc-work/consensusVariants/alignments/", genome, ".aligned.fasta", sep = "")

alignment = read.alignment(output, format = "FASTA", forceToLower = TRUE)

alignmentVector = vector(mode = "character", length = length(alignment$seq))

consensus = read.csv(paste("/rds/user/pw457/hpc-work/consensusVariants/whole_genome/cen180_consensus_", genome, ".csv", sep = ""))[1,2]
consensus = strsplit(consensus, split = "")[[1]]
consensus = tolower(consensus)

consensus.table = read.csv(paste("/rds/user/pw457/hpc-work/consensusVariants/whole_genome/cen180_frequencies_", genome, ".csv", sep = ""))

consensus.positions = which(consensus.table$in_consensus)

for(i in 1: length(alignmentVector))
{
  alignmentVector[i] = strsplit(alignment$seq[[i]], split = "")
}


chromosomes = unique(repeats$region.name)

files.done = list.files(path = ".", pattern = "variant.coordinates_*")

for(k in 1 : 5)
{
    if(length(grep(pattern = paste("variant.coordinates_", chromosomes[k], "_total", sep = ""), x =files.done)) < 1)
  {
    start.ID = min(which(repeats$region.name == chromosomes[k]))
    end.ID = max(which(repeats$region.name == chromosomes[k]))
    
    variants.total = 0
    coordinates = NULL
    
    for(i in start.ID : end.ID)
    {
      print(i)
      for(j in 1 : length(consensus.positions))
      {
        if(consensus[j] != alignmentVector[[i]][consensus.positions[j]])
        {
          variants.total = variants.total + 1
          ##save coordinate
          coordinates = c(coordinates, (repeats$start[i] + j))
        }
        alignmentVector[[i]][consensus.positions[j]] = "-"
      }
      for(j in 1 : length(consensus.positions))
      {
        if(alignmentVector[[i]][j] != "-")
        {
          variants.total = variants.total + 1
          ##save coordinate
          coordinates = c(coordinates, (repeats$start[i] + j))
        }
      }
    }
    write.csv(x = coordinates, file = paste("variant.coordinates_", chromosomes[k], ".csv", sep = ""), row.names = FALSE)
  } else
  {
    print("done before")
  }
}
























