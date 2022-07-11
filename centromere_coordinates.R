# load files with repeats
# check start sites, find outliers
# see how it plots


library(stringr)
library(base)
library(msa)
library(Biostrings)
library(seqinr)
library(kmer)

revCompString = function(DNAstr) {return(toupper(toString(reverseComplement(DNAString(DNAstr))))) }

columnRemove = function(DF, names){
  for(i in 1 : length(names))  {
    if(names[i] %in% names(DF)) {
      DF = DF[,-which(names[i] %in% names(DF))] 
    }
  }
  return(DF)
}

setwd("C:/Users/wlodz/Desktop/panC_analysis_athome/centromere_coordinates")



assemblies = list.files(path = "C:/Users/wlodz/Desktop/panC_analysis_athome/assemblies/66", pattern = ".f", full.names = T)
assemblies.names = list.files(path = "C:/Users/wlodz/Desktop/panC_analysis_athome/assemblies/66", pattern = ".f", full.names = F)


repeat.files = list.files(path = "C:/Users/wlodz/Desktop/panC_analysis_athome/truncated_repeats", pattern = "consensus.repetitiveness", full.names = TRUE)

repeat.files.accession = repeat.files
repeat.files.chr = repeat.files
for(i in 1 : length(repeat.files))
{  
  repeat.files.accession[i] = strsplit(repeat.files[i], split = "repetitiveness")[[1]][[2]]
  
  repeat.files.chr[i] = as.numeric(strsplit(repeat.files.accession[i], split = "")[[1]][length(strsplit(repeat.files.accession[i], split = "")[[1]]) - 4])
  
  temp = strsplit(repeat.files.accession[i], split = "[.]")[[1]]
  repeat.files.accession[i] = paste(temp[1 : (length(temp) - 2)], collapse = ".")
  
}

centromere.coordinates = data.frame(accession = vector(mode = "character", length = 330),
                                    chromosome.no = vector(mode = "character", length = 330),
                                    start = vector(mode = "character", length = 330),
                                    end = vector(mode = "character", length = 330))



ID = 0
repeats_per_array_min = 150
half_repeats_per_array_min = round(repeats_per_array_min/2, 0)
gaps_allowed = round(repeats_per_array_min/10,0) #extra repeats per array that can be gaps

for(i in 1 : length(assemblies))
{
  assembly.name = assemblies.names[i]
  assembly = read.fasta(file = assemblies[i])
  
  png(filename = paste(assembly.name, ".centromere.coordinates.png", sep = ""), width = 1000, height = 1000)
  par(mfrow = c(5,1))
  
  
  for(j in 1 : 5)
  {
    ID = ID + 1
    print(ID)
    chromosome = j
    repeats = read.csv(repeat.files[repeat.files.accession == assembly.name & repeat.files.chr == as.character(chromosome)])
    repeats$dist_to_prev = 0
    repeats$dist_to_next = 0
    repeats$out = TRUE
    
    for(k in 1 : nrow(repeats))
    {
      if(k - half_repeats_per_array_min > 0)
      {
        repeats$dist_to_prev[k] = repeats$start[k] - repeats$start[k - half_repeats_per_array_min]
      }
      if(k + half_repeats_per_array_min <= nrow(repeats))
      {
        repeats$dist_to_next[k] = repeats$start[k + half_repeats_per_array_min] - repeats$start[k]
      }
    }
    
    repeats$out[repeats$dist_to_next < (177*(half_repeats_per_array_min + gaps_allowed)) & 
                  repeats$dist_to_next > (177*(half_repeats_per_array_min - gaps_allowed))] = FALSE
    repeats$out[repeats$dist_to_prev < (177*(half_repeats_per_array_min + gaps_allowed)) & 
                  repeats$dist_to_prev > (177*(half_repeats_per_array_min - gaps_allowed))] = FALSE
    
    centromere.coordinates$accession[ID] = assembly.name
    centromere.coordinates$chromosome.no[ID] = j
    
    centromere.coordinates$start[ID] = min(repeats$start[!repeats$out]) - 250000
    centromere.coordinates$end[ID] = max(repeats$start[!repeats$out]) + 250000
    
    
    plot(x = repeats$start, y = seq(1, 1, length.out = length(repeats$start)), xlim = c(0,length(assembly[[j]])), ylim = c(0,2), pch = 20)
    points(x = repeats$start[!repeats$out], y = seq(1, 1, length.out = length(repeats$start[!repeats$out])), col = "red")
    abline(v = centromere.coordinates$start[ID])
    abline(v = centromere.coordinates$end[ID])
    
  }
  dev.off()
}
write.csv(x = centromere.coordinates, file = "2centromere.coordinates.csv")
