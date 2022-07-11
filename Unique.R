
library(stringr)
library(plot.matrix)
library("amap")
library("graph4lg")

names.of.repeats.files = list.files(path = "C:/Users/wlodz/Desktop/panC",pattern = "all.repeats.from.", recursive = TRUE, full.names = TRUE) 

cen180 = NA
for(i in 1 : length(names.of.repeats.files))
{
  print(i)
  repeats = read.csv(file = names.of.repeats.files[i])
  cen180 = rbind(cen180, repeats[!is.na(repeats$class),])
}

#repeats$chromosome = paste(repeats$repeats.genome, repeats$repeats.fasta.name, sep = "_")
#repeats = repeats[,-c(1,2,3)]
#repeats$width = repeats$repeats.end - repeats$repeats.start + 1
#names(repeats) = c("seq", "start", "end", "strand", "class", "chromosome", "width")
#cen180 = rbind(cen180, repeats)

cen180 = cen180[cen180$class == "aTha178",]
cen180$est.chromosome.no = ""

chr.names = which(str_detect(unique(cen180$chromosome), regex("Chr.")) | str_detect(unique(cen180$chromosome), regex("SUPER.")) | str_detect(unique(cen180$chromosome), regex("chr.")))
chr.names = unique(cen180$chromosome)[chr.names]
cen180b = cen180
cen180 = cen180[which(cen180$chromosome %in% chr.names),]

remove(repeats)

setwd("C:/Users/wlodz/Desktop/cenRepeats/unique/v4")

chr.1 = which(str_detect(unique(cen180$chromosome), regex("Chr1")) | str_detect(unique(cen180$chromosome), regex("SUPER_1")) | str_detect(unique(cen180$chromosome), regex("chr1")))
chr.1 = unique(cen180$chromosome)[chr.1]
cen180$est.chromosome.no[which(cen180$chromosome %in% chr.1)] = "chr1"
cen180.chr1 = cen180[which(cen180$chromosome %in% chr.1),]

chr.2 = which(str_detect(unique(cen180$chromosome), regex("Chr2")) | str_detect(unique(cen180$chromosome), regex("SUPER_2")) | str_detect(unique(cen180$chromosome), regex("chr2")))
chr.2 = unique(cen180$chromosome)[chr.2]
cen180$est.chromosome.no[which(cen180$chromosome %in% chr.2)] = "chr2"
cen180.chr2 = cen180[which(cen180$chromosome %in% chr.2),]

chr.3 = which(str_detect(unique(cen180$chromosome), regex("Chr3")) | str_detect(unique(cen180$chromosome), regex("SUPER_3")) | str_detect(unique(cen180$chromosome), regex("chr3")))
chr.3 = unique(cen180$chromosome)[chr.3]
cen180$est.chromosome.no[which(cen180$chromosome %in% chr.3)] = "chr3"
cen180.chr3 = cen180[which(cen180$chromosome %in% chr.3),]

chr.4 = which(str_detect(unique(cen180$chromosome), regex("Chr4")) | str_detect(unique(cen180$chromosome), regex("SUPER_4")) | str_detect(unique(cen180$chromosome), regex("chr4")))
chr.4 = unique(cen180$chromosome)[chr.4]
cen180$est.chromosome.no[which(cen180$chromosome %in% chr.4)] = "chr4"
cen180.chr4 = cen180[which(cen180$chromosome %in% chr.4),]

chr.5 = which(str_detect(unique(cen180$chromosome), regex("Chr5")) | str_detect(unique(cen180$chromosome), regex("SUPER_5")) | str_detect(unique(cen180$chromosome), regex("chr5")))
chr.5 = unique(cen180$chromosome)[chr.5]
cen180$est.chromosome.no[which(cen180$chromosome %in% chr.5)] = "chr5"
cen180.chr5 = cen180[which(cen180$chromosome %in% chr.5),]

write.csv(x = cen180, file = "allCEN180.panC.genomes.66.csv")
cen180 = read.csv(file = "allCEN180.panC.genomes.66.csv")

unique.matrix1 = matrix(ncol = length(unique(cen180.chr1$assembly)), nrow = length(unique(cen180.chr1$assembly)))
for(i in 1 : length(unique(cen180.chr1$assembly)))
{
  print(i)
  set1 = cen180.chr1$sequence_strand_adjusted[cen180.chr1$assembly == unique(cen180.chr1$assembly)[i]]
  for(ii in 1 : length(unique(cen180.chr1$assembly)))
  {
    set2 = cen180.chr1$sequence_strand_adjusted[cen180.chr1$assembly == unique(cen180.chr1$assembly)[ii]]
    unique.matrix1[i,ii] = length(which(set1 %in% set2)) /  (length(set1) + length(set2) / 2)
  }
}

tempcen180 = cen180.chr2
unique.matrixtemp = matrix(ncol = length(unique(tempcen180$assembly)), nrow = length(unique(tempcen180$assembly)))
for(i in 1 : length(unique(tempcen180$assembly)))
{
  print(i)
  set1 = tempcen180$sequence_strand_adjusted[tempcen180$assembly == unique(tempcen180$assembly)[i]]
  for(ii in 1 : length(unique(tempcen180$assembly)))
  {
    set2 = tempcen180$sequence_strand_adjusted[tempcen180$assembly == unique(tempcen180$assembly)[ii]]
    unique.matrixtemp[i,ii] = length(which(set1 %in% set2)) /  (length(set1) + length(set2) / 2)
  }
}
unique.matrix2 = unique.matrixtemp

tempcen180 = cen180.chr3
unique.matrixtemp = matrix(ncol = length(unique(tempcen180$assembly)), nrow = length(unique(tempcen180$assembly)))
for(i in 1 : length(unique(tempcen180$assembly)))
{
  print(i)
  set1 = tempcen180$sequence_strand_adjusted[tempcen180$assembly == unique(tempcen180$assembly)[i]]
  for(ii in 1 : length(unique(tempcen180$assembly)))
  {
    set2 = tempcen180$sequence_strand_adjusted[tempcen180$assembly == unique(tempcen180$assembly)[ii]]
    unique.matrixtemp[i,ii] = length(which(set1 %in% set2)) /  (length(set1) + length(set2) / 2)
  }
}
unique.matrix3 = unique.matrixtemp

tempcen180 = cen180.chr4
unique.matrixtemp = matrix(ncol = length(unique(tempcen180$assembly)), nrow = length(unique(tempcen180$assembly)))
for(i in 1 : length(unique(tempcen180$assembly)))
{
  print(i)
  set1 = tempcen180$sequence_strand_adjusted[tempcen180$assembly == unique(tempcen180$assembly)[i]]
  for(ii in 1 : length(unique(tempcen180$assembly)))
  {
    set2 = tempcen180$sequence_strand_adjusted[tempcen180$assembly == unique(tempcen180$assembly)[ii]]
    unique.matrixtemp[i,ii] = length(which(set1 %in% set2)) /  (length(set1) + length(set2) / 2)
  }
}
unique.matrix4 = unique.matrixtemp

tempcen180 = cen180.chr5
unique.matrixtemp = matrix(ncol = length(unique(tempcen180$assembly)), nrow = length(unique(tempcen180$assembly)))
for(i in 1 : length(unique(tempcen180$assembly)))
{
  print(i)
  set1 = tempcen180$sequence_strand_adjusted[tempcen180$assembly == unique(tempcen180$assembly)[i]]
  for(ii in 1 : length(unique(tempcen180$assembly)))
  {
    set2 = tempcen180$sequence_strand_adjusted[tempcen180$assembly == unique(tempcen180$assembly)[ii]]
    unique.matrixtemp[i,ii] = length(which(set1 %in% set2)) /  (length(set1) + length(set2) / 2)
  }
}
unique.matrix5 = unique.matrixtemp

write.table(unique.matrix1, "PanCnormalised.unique.matrix1.csv", row.names = unique(cen180.chr1$assembly), col.names = unique(cen180.chr1$assembly), sep = ",")
write.table(unique.matrix2, "PanCnormalised.unique.matrix2.csv", row.names = unique(cen180.chr2$assembly), col.names = unique(cen180.chr2$assembly), sep = ",")
write.table(unique.matrix3, "PanCnormalised.unique.matrix3.csv", row.names = unique(cen180.chr3$assembly), col.names = unique(cen180.chr3$assembly), sep = ",")
write.table(unique.matrix4, "PanCnormalised.unique.matrix4.csv", row.names = unique(cen180.chr4$assembly), col.names = unique(cen180.chr4$assembly), sep = ",")
write.table(unique.matrix5, "PanCnormalised.unique.matrix5.csv", row.names = unique(cen180.chr5$assembly), col.names = unique(cen180.chr5$assembly), sep = ",")

class(unique.matrix)

png("PanCnormalised.unique.matrix1.png", width = 1000, height = 1000, pointsize = 20)
plot(unique.matrix1, col = topo.colors, breaks = 100)
#plot(unique.matrix1, col=c('red', 'green'), breaks=c(1:100))
dev.off()

png("PanCnormalised.unique.matrix2.png", width = 1000, height = 1000, pointsize = 20)
plot(unique.matrix2, col = topo.colors, breaks = 100)
#plot(unique.matrix2, col=c('red', 'green'), breaks=c(1:100))
dev.off()

png("PanCnormalised.unique.matrix3.png", width = 1000, height = 1000, pointsize = 20)
plot(unique.matrix3, col = topo.colors, breaks = 100)
#plot(unique.matrix3, col=c('red', 'green'), breaks=c(1:100))
dev.off()

png("PanCnormalised.unique.matrix4.png", width = 1000, height = 1000, pointsize = 20)
plot(unique.matrix4, col = topo.colors, breaks = 100)
#plot(unique.matrix4, col=c('red', 'green'), breaks=c(1:100))
dev.off()

png("PanCnormalised.unique.matrix5.png", width = 1000, height = 1000, pointsize = 20)
plot(unique.matrix5, col = topo.colors, breaks = 100)
#plot(unique.matrix5, col=c('red', 'green'), breaks=c(1:100))
dev.off()



unique.matrix = matrix(ncol = length(unique(cen180$assembly)), nrow = length(unique(cen180$assembly)))
for(i in 1 : length(unique(cen180$assembly)))
{
  print(i)
  set1 = cen180$sequence_strand_adjusted[cen180$assembly == unique(cen180$assembly)[i]]
  for(ii in 1 : length(unique(cen180$assembly)))
  {
    set2 = cen180$sequence_strand_adjusted[cen180$assembly == unique(cen180$assembly)[ii]]
    unique.matrix[i,ii] = length(which(set1 %in% set2)) /  (length(set1) + length(set2) / 2)
  }
}
write.table(unique.matrix, "PanCnormalised.unique.matrix.csv", row.names = unique(cen180$assembly), col.names = unique(cen180$assembly), sep = ",")

png("PanCnormalised.unique.matrix.png", width = 1000, height = 1000, pointsize = 20)
plot(unique.matrix, col = topo.colors, breaks = 100)
dev.off()


colnames(unique.matrix) = 1:nrow(unique.matrix)
rownames(unique.matrix) = 1:nrow(unique.matrix)
symmetrix.Matrix = unique.matrix

for(i in 1 : nrow(symmetrix.Matrix))
{
  for(j in 1 : nrow(symmetrix.Matrix))
  {
    symmetrix.Matrix[i,j] = (symmetrix.Matrix[i,j] + symmetrix.Matrix[j,i]) / 2
    symmetrix.Matrix[j,i] = symmetrix.Matrix[i,j]
  }
}
symmetrix.Matrix = round(symmetrix.Matrix,3)

png("PanCsymmetric.unique.matrix.png", width = 1000, height = 1000, pointsize = 20)
plot(symmetrix.Matrix, col = topo.colors, breaks = 100)
dev.off()

symmetrix.Matrix[lower.tri(symmetrix.Matrix)] = t(symmetrix.Matrix)[lower.tri(symmetrix.Matrix)]
row.names(symmetrix.Matrix) = colnames(symmetrix.Matrix) = 1:nrow(unique.matrix)

cluster.Matrix = hcluster(symmetrix.Matrix)
reordered = reorder_mat(symmetrix.Matrix, as.character(cluster.Matrix$order))

orderWhole = cluster.Matrix$order

colnames(unique.matrix) = rownames(unique.matrix) = chr.names[as.integer(rownames(reordered))]

pdf("powTest_PanCreordered.symmetric.unique.matrix.pdf", width = 25, height = 25, pointsize = 25, onefile = TRUE)
for(i in 1 : 10)
{
  reorderedE = reordered^(i/10)
  plot(reorderedE, col = topo.colors, breaks = 100, main = paste("pow", i/10, " all chromosomes shared repeats to all repeats ratio", sep = ""), las = 2)
}
dev.off()

cluster.Matrix2 = hcluster(reordered)
pdf("PanCcluster dendrogram.pdf", width = 15, pointsize = 12)
plot(cluster.Matrix2, labels = unique(cen180$assembly)[as.integer(rownames(reordered))], main = "all chr dendrogram")
dev.off()

###########

symmetrix.Matrix = unique.matrix1
for(i in 1 : nrow(symmetrix.Matrix))
{
  for(j in 1 : nrow(symmetrix.Matrix))
  {
    symmetrix.Matrix[i,j] = (symmetrix.Matrix[i,j] + symmetrix.Matrix[j,i]) / 2
    symmetrix.Matrix[j,i] = symmetrix.Matrix[i,j]
  }
}
symmetrix.Matrix = round(symmetrix.Matrix,3)
symmetrix.Matrix[lower.tri(symmetrix.Matrix)] = t(symmetrix.Matrix)[lower.tri(symmetrix.Matrix)]
row.names(symmetrix.Matrix) = colnames(symmetrix.Matrix) = 1:nrow(symmetrix.Matrix)
cluster.Matrix = hcluster(symmetrix.Matrix)
reordered = reorder_mat(symmetrix.Matrix, as.character(cluster.Matrix$order))
chr.single.names = chr.names[(str_detect(chr.names, regex(".hr1")))]
pdf("PanCreordered.symmetric.unique.matrix.1.pdf", width = 25, height = 25, pointsize = 25)
plot(reordered, col = topo.colors, breaks = 100, main = "chr1 shared repeats to all repeats ratio", las = 2)
dev.off()


cluster.Matrix2 = hcluster(reordered)
colnames(reordered) = rownames(reordered) = unique(cen180$assembly)[as.integer(rownames(reordered))]
pdf("PanCcluster dendrogram.1.pdf", width = 15, pointsize = 12)
plot(cluster.Matrix2, labels = rownames(reordered), main = "chr1 dendrogram")
dev.off()
######################

symmetrix.Matrix = unique.matrix2
for(i in 1 : nrow(symmetrix.Matrix))
{
  for(j in 1 : nrow(symmetrix.Matrix))
  {
    symmetrix.Matrix[i,j] = (symmetrix.Matrix[i,j] + symmetrix.Matrix[j,i]) / 2
    symmetrix.Matrix[j,i] = symmetrix.Matrix[i,j]
  }
}
symmetrix.Matrix = round(symmetrix.Matrix,3)
symmetrix.Matrix[lower.tri(symmetrix.Matrix)] = t(symmetrix.Matrix)[lower.tri(symmetrix.Matrix)]
row.names(symmetrix.Matrix) = colnames(symmetrix.Matrix) = 1:nrow(symmetrix.Matrix)
cluster.Matrix = hcluster(symmetrix.Matrix)
reordered = reorder_mat(symmetrix.Matrix, as.character(cluster.Matrix$order))
chr.single.names = chr.names[(str_detect(chr.names, regex(".hr2")))]
pdf("PanCreordered.symmetric.unique.matrix.2.pdf", width = 25, height = 25, pointsize = 25)
plot(reordered, col = topo.colors, breaks = 100, main = "chr2 shared repeats to all repeats ratio", las = 2)
dev.off()


cluster.Matrix2 = hcluster(reordered)
colnames(reordered) = rownames(reordered) = unique(cen180$assembly)[as.integer(rownames(reordered))]
pdf("PanCcluster dendrogram.2.pdf", width = 15, pointsize = 12)
plot(cluster.Matrix2, labels = rownames(reordered), main = "chr2 dendrogram")
dev.off()

######################

symmetrix.Matrix = unique.matrix3
for(i in 1 : nrow(symmetrix.Matrix))
{
  for(j in 1 : nrow(symmetrix.Matrix))
  {
    symmetrix.Matrix[i,j] = (symmetrix.Matrix[i,j] + symmetrix.Matrix[j,i]) / 2
    symmetrix.Matrix[j,i] = symmetrix.Matrix[i,j]
  }
}
symmetrix.Matrix = round(symmetrix.Matrix,3)
symmetrix.Matrix[lower.tri(symmetrix.Matrix)] = t(symmetrix.Matrix)[lower.tri(symmetrix.Matrix)]
row.names(symmetrix.Matrix) = colnames(symmetrix.Matrix) = 1:nrow(symmetrix.Matrix)
cluster.Matrix = hcluster(symmetrix.Matrix)
reordered = reorder_mat(symmetrix.Matrix, as.character(cluster.Matrix$order))
chr.single.names = chr.names[(str_detect(chr.names, regex(".hr3")))]
pdf("PanCreordered.symmetric.unique.matrix.3.pdf", width = 25, height = 25, pointsize = 25)
plot(reordered, col = topo.colors, breaks = 100, main = "chr3 shared repeats to all repeats ratio", las = 2)
dev.off()


cluster.Matrix2 = hcluster(reordered)
colnames(reordered) = rownames(reordered) = unique(cen180$assembly)[as.integer(rownames(reordered))]
pdf("PanCcluster dendrogram.3.pdf", width = 15, pointsize = 12)
plot(cluster.Matrix2, labels = rownames(reordered), main = "chr3 dendrogram")
dev.off()

######################

symmetrix.Matrix = unique.matrix4
for(i in 1 : nrow(symmetrix.Matrix))
{
  for(j in 1 : nrow(symmetrix.Matrix))
  {
    symmetrix.Matrix[i,j] = (symmetrix.Matrix[i,j] + symmetrix.Matrix[j,i]) / 2
    symmetrix.Matrix[j,i] = symmetrix.Matrix[i,j]
  }
}
symmetrix.Matrix = round(symmetrix.Matrix,3)
symmetrix.Matrix[lower.tri(symmetrix.Matrix)] = t(symmetrix.Matrix)[lower.tri(symmetrix.Matrix)]
row.names(symmetrix.Matrix) = colnames(symmetrix.Matrix) = 1:nrow(symmetrix.Matrix)
cluster.Matrix = hcluster(symmetrix.Matrix)
reordered = reorder_mat(symmetrix.Matrix, as.character(cluster.Matrix$order))
chr.single.names = chr.names[(str_detect(chr.names, regex(".hr4")))]
pdf("PanCreordered.symmetric.unique.matrix.4.pdf", width = 25, height = 25, pointsize = 25)
plot(reordered, col = topo.colors, breaks = 100, main = "chr4 shared repeats to all repeats ratio", las = 2)
dev.off()


cluster.Matrix2 = hcluster(reordered)
colnames(reordered) = rownames(reordered) = unique(cen180$assembly)[as.integer(rownames(reordered))]
pdf("PanCcluster dendrogram.4.pdf", width = 15, pointsize = 12)
plot(cluster.Matrix2, labels = rownames(reordered), main = "chr4 dendrogram")
dev.off()

######################

symmetrix.Matrix = unique.matrix5
for(i in 1 : nrow(symmetrix.Matrix))
{
  for(j in 1 : nrow(symmetrix.Matrix))
  {
    symmetrix.Matrix[i,j] = (symmetrix.Matrix[i,j] + symmetrix.Matrix[j,i]) / 2
    symmetrix.Matrix[j,i] = symmetrix.Matrix[i,j]
  }
}
symmetrix.Matrix = round(symmetrix.Matrix,3)
symmetrix.Matrix[lower.tri(symmetrix.Matrix)] = t(symmetrix.Matrix)[lower.tri(symmetrix.Matrix)]
row.names(symmetrix.Matrix) = colnames(symmetrix.Matrix) = 1:nrow(symmetrix.Matrix)
cluster.Matrix = hcluster(symmetrix.Matrix)
reordered = reorder_mat(symmetrix.Matrix, as.character(cluster.Matrix$order))
chr.single.names = chr.names[(str_detect(chr.names, regex(".hr5")))]
pdf("PanCreordered.symmetric.unique.matrix.5.pdf", width = 25, height = 25, pointsize = 25)
plot(reordered, col = topo.colors, breaks = 100, main = "chr5 shared repeats to all repeats ratio", las = 2)
dev.off()


cluster.Matrix2 = hcluster(reordered)
colnames(reordered) = rownames(reordered) = unique(cen180$assembly)[as.integer(rownames(reordered))]
pdf("PanCcluster dendrogram.5.pdf", width = 15, pointsize = 12)
plot(cluster.Matrix2, labels = rownames(reordered), main = "chr5 dendrogram")
dev.off()


##########################
#
#
#   all chromosomes separately, ordered by chromosome
#
#
#########################

tempcen180 = cen180


symmetrix.Matrix[lower.tri(symmetrix.Matrix)] = t(symmetrix.Matrix)[lower.tri(symmetrix.Matrix)]
row.names(symmetrix.Matrix) = colnames(symmetrix.Matrix) = 1:nrow(unique.matrix)
cluster.Matrix = hcluster(symmetrix.Matrix)
reordered = reorder_mat(symmetrix.Matrix, as.character(cluster.Matrix$order))
colnames(unique.matrix) = rownames(unique.matrix) = unique(cen180$assembly)[as.integer(rownames(reordered))]
#reorder tempcen180 by assembly accordinf to the orderWhole
tempcen180 = NULL
for(i in 1 : length(unique(cen180$assembly)))
{
  print(i)
  tempcen180 = rbind(tempcen180, cen180[cen180$assembly == unique(cen180$assembly)[orderWhole[i]],])
}
#reorder tempcen180 by chromosome 1:5
tempcen180b = tempcen180
tempcen180 = NULL
for(i in 1 : 5)
{
  print(i)
  tempcen180 = rbind(tempcen180, tempcen180b[tempcen180b$est.chromosome.no == paste("chr", i, sep = ""),])
}

unique.matrixtemp = matrix(ncol = (length(unique(tempcen180$assembly)) * 5), nrow = (length(unique(tempcen180$assembly)) * 5))

list = data.frame(genome = vector(length = nrow(unique.matrixtemp), mode = "character"), chr = vector(length = nrow(unique.matrixtemp), mode = "character"), merge = vector(length = nrow(unique.matrixtemp), mode = "character"))

for(i in 1 : nrow(unique.matrixtemp))
{
  print(i)
  if(i%%66 == 0)
  {
    list$genome[i] = unique(tempcen180$assembly)[66]
  } else
  {
    list$genome[i] = unique(tempcen180$assembly)[i%%66]
  }
  list$chr[i] = paste("chr", (((i-1)%/%66)+1), sep = "")
  list$merge[i] = paste(list$genome[i], list$chr[i], sep = "_")
}

for(i in 1 : nrow(unique.matrixtemp))
{
  print(i)
  set1 = tempcen180$sequence_strand_adjusted[tempcen180$assembly == list$genome[i] & tempcen180$est.chromosome.no == list$chr[i]]
  for(ii in 1 : nrow(unique.matrixtemp))
  {
    set2 = tempcen180$sequence_strand_adjusted[tempcen180$assembly == list$genome[ii] & tempcen180$est.chromosome.no == list$chr[ii]]
    unique.matrixtemp[i,ii] = length(which(set1 %in% set2)) / (length(set1) + length(set2) / 2)
  }
}

symmetrix.Matrix = unique.matrixtemp
for(i in 1 : nrow(symmetrix.Matrix))
{
  for(j in 1 : nrow(symmetrix.Matrix))
  {
    symmetrix.Matrix[i,j] = (symmetrix.Matrix[i,j] + symmetrix.Matrix[j,i]) / 2
    symmetrix.Matrix[j,i] = symmetrix.Matrix[i,j]
  }
}
symmetrix.Matrix = round(symmetrix.Matrix,3)
symmetrix.Matrix[lower.tri(symmetrix.Matrix)] = t(symmetrix.Matrix)[lower.tri(symmetrix.Matrix)]
row.names(symmetrix.Matrix) = colnames(symmetrix.Matrix) = 1:nrow(symmetrix.Matrix)
#cluster.Matrix = hcluster(symmetrix.Matrix)
#reordered = reorder_mat(symmetrix.Matrix, as.character(cluster.Matrix$order))
pdf("PanCreordered.symmetruc.unique.matrix.all.pdf", width = 100, height = 100, pointsize = 25)
plot(symmetrix.Matrix, col = topo.colors, breaks = 100, main = "all chromosomes shared repeats to all repeats ratio", las = 2)
dev.off()

##
#code below doesnt work properly#
##
cluster.Matrix2 = hcluster(symmetrix.Matrix)
cluster.Matrix2 = reorder(cluster.Matrix2, wts = cluster.Matrix2$order)
for(i in 1 : length(cluster.Matrix2$order))
colnames(symmetrix.Matrix) = rownames(symmetrix.Matrix) = list$merge
pdf("PanCcluster dendrogram.all.pdf", width = 60, height = 28, pointsize = 12)
plot(cluster.Matrix2, labels = list$merge, main = "all chromosomes dendrogram")
plot(cluster.Matrix2, labels = FALSE)
dev.off()

write.csv(x = list$merge, file = "list.csv")

#######################################################################
#         how many repeats from chromosome A (set1) are present       #
#        on chromosome B (set2) using exact matching with %in%        #
#######################################################################
setwd("C:/Users/wlodz/Desktop/IansUnique/") #set working directory where the file with repeats is ("all.repeats.weigel.sanger.colcen.140821.csv")
all.repeats = read.csv("all.repeats.weigel.sanger.colcen.140821.csv")
all.cen180 = all.repeats[all.repeats$repeats.class == "CEN180",]
remove(all.repeats)

all.cen180$merge.name = paste(all.cen180$repeats.genome, all.cen180$repeats.fasta.name, sep = "_")
all.cen180$is.chromosome = FALSE
all.cen180$is.chromosome = TRUE * (substr(all.cen180$repeats.fasta.name,1,3) == "Chr")
all.cen180.chr = all.cen180[all.cen180$is.chromosome == TRUE,]
file.names2 = unique(all.cen180.chr$merge.name)#get a list of all chromosomes to pairwise compare
remove(all.cen180.chr)

unique.matrix.chr = matrix(ncol = length(file.names2), nrow = length(file.names2)) #prepare the matrix to put the values in

for(i in 1 : length(file.names2)) #iterate over all unique chromosomes
{
  print(i)
  set1 = all.cen180$repeats.sequence[all.cen180$merge.name == file.names2[i]] #prepare set1 that will be compared with...
  for(ii in 1 : length(file.names2))  # ...all other chromosomes...
  {
    set2 = all.cen180$repeats.sequence[all.cen180$merge.name == file.names2[ii]] #... as a set2
    unique.matrix.chr[i,ii] = length(which(set1 %in% set2))#... save number of matched repeats from set1 in set2
  }
}

write.table(unique.matrix.chr, "unique.matrix.chr.csv", row.names = file.names2, col.names = file.names2, sep = ",")
###########################
rownames(unique.matrix.chr) = file.names2
colnames(unique.matrix.chr) = file.names2

rbPal <- colorRampPalette(c('red','green'))

png("unique.matrix.chr.png", width = 2000, height = 2000, pointsize = 40)
plot(unique.matrix.chr, col = rbPal(4), breaks = c(0,100,1000,15000))
abline(h = (seq(0, length(file.names), 5) + 0.5), v = (seq(0, length(file.names), 5) + 0.5), col = "blue", lwd = 5)
dev.off()


unique.matrix.chr.log = log10(unique.matrix.chr)
unique.matrix.chr.log[is.infinite(unique.matrix.chr.log)] = 0
hist(unique.matrix.chr.log)

breaks = seq(0, 5, 0.2)
png("unique.matrix.chr.log.png", width = 2000, height = 2000, pointsize = 40)
plot(unique.matrix.chr.log, col = rbPal(length(breaks)), breaks = breaks)
abline(h = (seq(0, length(file.names), 5) + 0.5), v = (seq(0, length(file.names), 5) + 0.5), col = "blue", lwd = 5)
#abline(h = (seq(0, length(file.names), 5)), v = (seq(0, length(file.names), 5)))
dev.off()





list = 0
for(i in 1 : 115)
{
  list[i] = strsplit(file.names2, "_")[[i]][2]
}

unique.chr1 = unique.matrix.chr[list == "Chr1",list == "Chr1"]
unique.chr2 = unique.matrix.chr[list == "Chr2",list == "Chr2"]
unique.chr3 = unique.matrix.chr[list == "Chr3",list == "Chr3"]
unique.chr4 = unique.matrix.chr[list == "Chr4",list == "Chr4"]
unique.chr5 = unique.matrix.chr[list == "Chr5",list == "Chr5"]

unique.genomes = unique.chr1 + unique.chr2 + unique.chr3 + unique.chr4 + unique.chr5

unique.genomes.log = log10(unique.genomes)
breaks = seq(0, 5, 0.2)
png("unique.genomes.log.png", width = 1000, height = 1000, pointsize = 20)
plot(unique.genomes.log, col = rbPal(length(breaks)), breaks = breaks)
#abline(h = (seq(0, length(file.names), 5) + 0.5), v = (seq(0, length(file.names), 5) + 0.5), col = "blue", lwd = 5)
#abline(h = (seq(0, length(file.names), 5)), v = (seq(0, length(file.names), 5)))
dev.off()

unique.genomes.norm = unique.genomes
for(i in 1 : 23)
{
  for(ii in 1 : 23)
  {
    if(i != ii)
    {
      unique.genomes.norm[i,ii] = unique.genomes.norm[i,ii] / unique.genomes.norm[i,i]
    }
  }
  unique.genomes.norm[i,i] = 1
}
unique.genomes.norm = unique.genomes.norm * 100

for(i in 1 : 22)
{
  for(ii in (i + 1) : 23)
  {
    unique.genomes.norm[i,ii] = unique.genomes.norm[i,ii] + unique.genomes.norm[ii,i]
    unique.genomes.norm[ii,i] = 0
  }
}

unique.genomes.norm[unique.genomes.norm == 0] = NA

breaks = seq(0, 100, 5)

png("unique.genomes.norm.png")
plot(unique.genomes.norm, col = rbPal(length(breaks)), breaks = breaks, na.cell = FALSE, xaxt='n', yaxt='n', ann=FALSE)
#abline(h = (seq(0, length(file.names), 5) + 0.5), v = (seq(0, length(file.names), 5) + 0.5), col = "blue", lwd = 5)
for(i in 1 : 23)
{
  text(x = i+1, y = 24-i, colnames(unique.genomes.norm)[i])
}
dev.off()










