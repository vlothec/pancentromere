library(stringr)
library(seqinr)
library(base)

setwd("C:/Users/wlodz/Desktop/cenRepeats/unique")

repeats.files = list.files(path = ".",pattern = "all.repeats.from.", recursive = TRUE) 

divide.and.fill.names = function(row, uniques)
{
  name.temp = unlist(strsplit(row[7], split = "_"))
  return(list(name.temp[length(name.temp)], which(uniques == row[7])))
}

for(i in 1 : length(repeats.files))
{
  print(paste("Working on file ", i, ": ", repeats.files[i], sep = ""))
  
  repeats.temp = read.csv(repeats.files[i], header = TRUE)
  
  repeats.temp = cbind(repeats.temp, data.frame(t(sapply(apply(repeats.temp, 1, divide.and.fill.names, unique(repeats.temp$region.name)),c))))
  
  repeats.temp$X1 = unlist(repeats.temp$X1)
  repeats.temp$X2 = unlist(repeats.temp$X2)
  
  names(repeats.temp)[8] = "chromosome"
  names(repeats.temp)[9] = "est.chromosome.no"
  
  name.temp = strsplit(repeats.temp$region.name[1], split = "_")[[1]]
  repeats.temp$fasta.file.name = paste(name.temp[1 : (length(name.temp) - 1)], collapse = "_")
  remove(name.temp)
  
  repeats.temp = repeats.temp[, -which(names(repeats.temp) == "region.name")]
  
  write.csv(x = repeats.temp, file = paste("ed11.01.22.", repeats.files[i], sep = ""), row.names = FALSE)
  
  remove(repeats.temp)
  gc()
}
remove(repeats.files, i, divide.and.fill.names)
