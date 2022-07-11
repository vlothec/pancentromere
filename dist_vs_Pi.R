library(stringr)
library(base)
library(msa)
library(Biostrings)
library(seqinr)
library(microseq)
library(R.utils)
library(data.table)
library(plot.matrix)



setwd("C:/Users/wlodz/Desktop/temp/dist_vs_PI")
###################################### Geo distance vs PI

PI = read.table(file = "Pi_Chr1to5_66Athaliana.txt", header = TRUE)
distances = read.csv(file = "distances_66acc.csv", header = TRUE, row.names = 1)
summary.table = read.csv(file = "summary_table_accessions.csv")
names(distances) %in%  unique(PI$acc1)

for(i in 1 : length(unique(PI$acc1)))
{
  single = PI$acc1[PI$acc1 == unique(PI$acc1)[i]][1]
  if(length(grep(single, summary.table$accession_id)) == 1)
  {
    PI$acc1[PI$acc1 == unique(PI$acc1)[i]] = summary.table$accession_name[grep(single, summary.table$accession_id)]
  }
}
for(i in 1 : length(unique(PI$acc2)))
{
  single = PI$acc2[PI$acc2 == unique(PI$acc2)[i]][1]
  if(length(grep(single, summary.table$accession_id)) == 1)
  {
    PI$acc2[PI$acc2 == unique(PI$acc2)[i]] = summary.table$accession_name[grep(single, summary.table$accession_id)]
  }
}
PI$acc1[PI$acc1 == "Kew"] = "ddAraThal4"
PI$acc2[PI$acc2 == "Kew"] = "ddAraThal4"
PI$acc1[PI$acc1 == "11C1"] = "X11C1"
PI$acc2[PI$acc2 == "11C1"] = "X11C1"



unique(PI$acc1)
distances = as.matrix(distances)
for(i in 1 : length(colnames(distances)))
{
  colnames(distances)[i] = paste(strsplit(colnames(distances)[i], split = "[.]")[[1]], collapse = "_")
  rownames(distances)[i] = paste(strsplit(rownames(distances)[i], split = "[.]")[[1]], collapse = "_")
}
for(i in 1 : length(unique(PI$acc1)))
{
  PI$acc1[PI$acc1 == unique(PI$acc1)[i]] = paste(strsplit(unique(PI$acc1)[i], split = "[-]")[[1]], collapse = "_")
  PI$acc2[PI$acc2 == unique(PI$acc2)[i]] = paste(strsplit(unique(PI$acc2)[i], split = "[-]")[[1]], collapse = "_")
}


PI$acc1[PI$acc1 == "Ler_0"] = "Ler_0_110x"
PI$acc2[PI$acc2 == "Ler_0"] = "Ler_0_110x"

unique(PI$acc1)[!(unique(PI$acc1) %in% colnames(distances))]


PI$geo = 0
for(ID in 1 : nrow(PI))
{
  print(ID)
  PI$geo[ID] = distances[which(PI$acc1[ID] == colnames(distances)),which(PI$acc2[ID] == colnames(distances))]
}

distances = distances[,order(colnames(distances))]
distances = distances[order(rownames(distances)),]
write.csv(x = PI, file = "PI.with.geo.csv")

pdf(file = "geo vs pi.pdf")
plot(PI$geo, PI$pi, pch = 20)
dev.off()

