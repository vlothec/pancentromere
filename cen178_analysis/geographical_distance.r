
.libPaths("/home/pw457/Rpackages")

library(stringr)
library(base)
library(seqinr)#
library(ggplot2)
library("geosphere")



setwd("/rds/user/pw457/hpc-work/geo_distance")
accessions = read.csv(file = "summary_table_accessions.csv")

accessions = accessions[accessions$Column1]

distances = matrix(data = 0,nrow = nrow(accessions), ncol = nrow(accessions))
colnames(distances) = accessions$accession_name
rownames(distances) = accessions$accession_name

for(i in 1 : nrow(distances))
{
  print(i)
  for(j in 1 : nrow(distances))
  {
    p1 = c(accessions$longitude[i], accessions$latitude[i])
    p2 = c(accessions$longitude[j], accessions$latitude[j])
    
    distances[i,j] = distGeo(p1, p2, a=6378137, f=1/298.257223563)/1000
  }
}

write.csv(distances, file = "distances_66acc.csv")


