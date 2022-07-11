library(stringr)
library(base)
library(msa)
library(Biostrings)
library(seqinr)
library(Matrix)

#Figure 1 - plot of pairwise CEN178 sharing vs geographic distance

setwd("C:/Users/wlodz/Desktop/panC_analysis_athome/sharing_vs_geo")


#fix summary table
summary.table = read.csv("summary_table_accessions.csv")
for(i in 1 : nrow(summary.table))
{
  summary.table$accession_id[i] = strsplit(summary.table$accession_id[i], split = "[,]")[[1]][1]
  summary.table$accession_id[i] = paste(strsplit(summary.table$accession_id[i], split = "[-]")[[1]], collapse = ".")
}


#fix distance values names

distance.values = read.csv(file = "distances_66acc.csv", row.names = 1)
for(i in 1 : length(names(distance.values)))
{
  names(distance.values)[i] = paste(strsplit(names(distance.values)[i], split = "[.]")[[1]], collapse = "-")
}

names(distance.values)[5] = "11C1"

summary.table$accession_name %in% names(distance.values)

#fix sharing values names

sharing.values = read.csv(file = "PanCnormalised.unique.matrix.csv") ##TODO change it to a new one
for(i in 1 : length(names(sharing.values)))
{
  if(length(strsplit(names(sharing.values)[i], split = "X")[[1]]) > 1)
  {
    names(sharing.values)[i] = strsplit(names(sharing.values)[i], split = "X")[[1]][2]
  }
}
names(sharing.values)[39] = "ddAraThal4.1.primary.fa"
summary.table$accession_id %in% names(sharing.values)

for(i in 1 : length(names(sharing.values)))
{
  names(sharing.values)[i] = summary.table$accession_name[summary.table$accession_id == names(sharing.values)[i]]
}
summary.table$accession_name %in% names(sharing.values)


names(distance.values) %in% names(sharing.values)

#reorder distances
rownames(distance.values) = names(distance.values)
colnames(distance.values) = names(distance.values)
names(distance.values)

distance.values = distance.values[summary.table$accession_name , summary.table$accession_name]

#reorder sharing
rownames(sharing.values) = names(sharing.values)
colnames(sharing.values) = names(sharing.values)
names(sharing.values)
sharing.values = sharing.values[summary.table$accession_name , summary.table$accession_name]

####make the table

sharing_vs_geo = data.frame(nameA = NULL, nameB = NULL, sharing = NULL, geo = NULL)
for(i in 1 : 65)
{
  print(i)
  for(j in (i + 1) : 66)
  {
    nameA = summary.table$accession_name[i]
    nameB = summary.table$accession_name[j]
    sharing = sharing.values[i,j]
    geo = distance.values[i,j]
    sharing_vs_geo = rbind(sharing_vs_geo, data.frame(nameA, nameB, sharing, geo))
  }
}
sharing_vs_geo$sharing = (sharing_vs_geo$sharing/2)*3


ggplot(data = sharing_vs_geo, aes(sharing,geo)) + 
  geom_bin2d(bins = 100) + 
  theme_bw() + 
  scale_fill_gradient(low = "black", high = "red") +
  xlab("shared repeats/average number") +
  ylab("geographical distance, km") + 
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text.x = element_text(size = 14),  axis.text.y = element_text(size = 14))

ggsave(filename = "bin2d sharing vs geo.pdf", width = 8, height = 7)


