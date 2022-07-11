library(stringr)
library(base)
library(msa)
library(Biostrings)
library(seqinr)
library(kmer)


setwd("C:/Users/wlodz/Desktop/panC_analysis_athome/detect_large_duplications")

repeat.files = list.files(path = "C:/Users/wlodz/Desktop/panC_analysis_athome/truncated_repeats", pattern = "consensus.repetitiveness", full.names = TRUE)

HOR.files = list.files(path = "C:/Users/wlodz/Desktop/panC_analysis_athome/detect_large_duplications/t0c5", pattern = "HORs_")

repeat.files.accession = repeat.files
HOR.files.accession = HOR.files
HOR.files.chr = HOR.files
repeat.files.chr = repeat.files

for(i in 1 : length(HOR.files))
{
  
  HOR.files.accession[i] = strsplit(HOR.files[i], split = "from.")[[1]][[2]]
  HOR.files.accession[i] = paste(strsplit(HOR.files.accession[i], split = "")[[1]][1:(length(strsplit(HOR.files.accession[i], split = "")[[1]])-8)], collapse = "")
  if(HOR.files.accession[i] == "ddAraThal4.1.primary.fa")
  {
    HOR.files.chr[i] = as.numeric(strsplit(HOR.files[i], split = "_")[[1]][[3]])
  } else
  {
    temp = strsplit(HOR.files[i], split = "_")[[1]][[2]]
    HOR.files.chr[i] = as.numeric(strsplit(temp, split = "r")[[1]][[2]])
  }
}

HOR.files = list.files(path = "C:/Users/wlodz/Desktop/panC_analysis_athome/detect_large_duplications/t0c5", pattern = "HORs_", full.names = TRUE)

for(i in 1 : length(repeat.files))
{  
  repeat.files.accession[i] = strsplit(repeat.files[i], split = "repetitiveness")[[1]][[2]]
  
  repeat.files.chr[i] = as.numeric(strsplit(repeat.files.accession[i], split = "")[[1]][length(strsplit(repeat.files.accession[i], split = "")[[1]]) - 4])
  
  temp = strsplit(repeat.files.accession[i], split = "[.]")[[1]]
  repeat.files.accession[i] = paste(temp[1 : (length(temp) - 2)], collapse = ".")
  
}


HOR.files.accession %in% repeat.files.accession

distances = vector(mode = "list", length = length(HOR.files))
  
for(i in 1 : length(HOR.files))
{
  chromosome = HOR.files.chr[i]
  assembly = HOR.files.accession[i]
  print(paste(i, ": working on chromosome ", chromosome, " from ", assembly, sep = ""))
  hors = read.csv(HOR.files[i])
  #repeats = read.csv(file = repeat.files[which(repeat.files.accession == assembly & repeat.files.chr == chromosome)])
  
  hors = hors[hors$total_variant <= 0,]
  hors = hors[hors$block.size.in.units >= 5,]
  #hors = hors[abs(hors$start_A - hors$start_B) >= 400,]
  
  distances[i] = list(abs(hors$start_A - hors$start_B))
}




pdf(file = "boxplot distances.pdf", width = 30, height = 20)
par(mfrow = c(1,5))

stripchart(distances[HOR.files.chr == "1"], method="jitter", pch = 20, xlim = c(0,20000))
stripchart(distances[HOR.files.chr == "2"], method="jitter", pch = 20, xlim = c(0,20000))
stripchart(distances[HOR.files.chr == "3"], method="jitter", pch = 20, xlim = c(0,20000))
stripchart(distances[HOR.files.chr == "4"], method="jitter", pch = 20, xlim = c(0,20000))
stripchart(distances[HOR.files.chr == "5"], method="jitter", pch = 20, xlim = c(0,20000))

dev.off()

#

#read a random HOR file fisrt to iniate this table   #######
hors$duplication_id = 0
hors$accession = ""

duplicated.hors = hors[0,]

duplication_ID = 1

for(i in 1 : length(HOR.files))
{
  chromosome = HOR.files.chr[i]
  assembly = HOR.files.accession[i]
  print(paste(i, ": working on chromosome ", chromosome, " from ", assembly, sep = ""))
  hors = read.csv(HOR.files[i])
  #repeats = read.csv(file = repeat.files[which(repeat.files.accession == assembly & repeat.files.chr == chromosome)])
  
  hors = hors[!(hors$start.A.bp  %in%  boxplot.stats(hors$start.A.bp)$out),]
  hors = hors[!(hors$start.B.bp  %in%  boxplot.stats(hors$start.B.bp)$out),]
  
  hors = hors[hors$total_variant <= 0,]
  hors = hors[hors$block.size.in.units >= 5,]
  
  
  hors = hors[abs(hors$start_A - hors$start_B) >= 561,] #561 repeats is more or less 100 kb
  
  
  if(nrow(hors) > 0)
  {
    hors$accession = assembly
    hist = hist(abs(hors$start_A - hors$start_B), breaks = seq(0,22000, 50), plot = FALSE)
    
    breaks = NULL
    broken = FALSE
    
    for(j in 1 : length(hist$counts))
    {
      if(hist$counts[j] == 0)
      {
        if(broken)
        {
          
        } else
        {
          breaks = c(breaks, j)
          broken = TRUE
        }
      } else
      {
        broken = FALSE
      }
    }
    
    
    for(j in 1 : length(breaks))
    {
      if(j == length(breaks))
      {
        duplicated.hors.t = hors[abs(hors$start_A - hors$start_B) >= hist$breaks[breaks[j]],]
        if(nrow(duplicated.hors.t) > 0)
        {
          duplicated.hors.t = duplicated.hors.t[!(duplicated.hors.t$start.A.bp %in% boxplot.stats(duplicated.hors.t$start.A.bp)$out),]
          duplicated.hors.t$duplication_id = duplication_ID
          duplication_ID = duplication_ID + 1
          duplicated.hors = rbind(duplicated.hors, duplicated.hors.t)
        }
      } else
      {
        duplicated.hors.t = hors[abs(hors$start_A - hors$start_B) >= hist$breaks[breaks[j]] & abs(hors$start_A - hors$start_B) < hist$breaks[breaks[j+1]],]
        if(nrow(duplicated.hors.t) > 0)
        {
          duplicated.hors.t = duplicated.hors.t[!(duplicated.hors.t$start.A.bp %in% boxplot.stats(duplicated.hors.t$start.A.bp)$out),]
          duplicated.hors.t$duplication_id = duplication_ID
          duplication_ID = duplication_ID + 1
          duplicated.hors = rbind(duplicated.hors, duplicated.hors.t)
        }
      }
    }
  }
}

all.hor.suplications = duplicated.hors
write.csv(x = all.hor.suplications, file = "all.hor.duplications.csv")

### ### ### ### ### 
### 
### manual edits based on the plots
### 
if(FALSE)
{
  all.hor.suplications = read.csv(file = "all.hor.duplications.csv")
  
  #split these
  
  #24
  all.HORdups.temp = all.hor.suplications[all.hor.suplications$duplication_id == 24,]
  all.hor.suplications$duplication_id[all.hor.suplications$duplication_id == 24 & all.hor.suplications$start.A.bp >= 15102157] = max(all.hor.suplications$duplication_id) + 1
  
  #438
  all.HORdups.temp = all.hor.suplications[all.hor.suplications$duplication_id == 438,]
  all.hor.suplications$duplication_id[all.hor.suplications$duplication_id == 438 & all.hor.suplications$start.A.bp >= 7131852] = max(all.hor.suplications$duplication_id) + 1
  
  #624
  all.HORdups.temp = all.hor.suplications[all.hor.suplications$duplication_id == 624,]
  all.hor.suplications$duplication_id[all.hor.suplications$duplication_id == 624 & all.hor.suplications$start.A.bp >= 15991645] = max(all.hor.suplications$duplication_id) + 1
  
  #162
  all.HORdups.temp = all.hor.suplications[all.hor.suplications$duplication_id == 162,]
  all.hor.suplications$duplication_id[all.hor.suplications$duplication_id == 162 & all.hor.suplications$start.A.bp >= 17435159] = max(all.hor.suplications$duplication_id) + 1
  
  all.HORdups.temp = all.hor.suplications[all.hor.suplications$duplication_id == max(all.hor.suplications$duplication_id),]
  all.hor.suplications$duplication_id[all.hor.suplications$duplication_id == max(all.hor.suplications$duplication_id) & all.hor.suplications$start.A.bp >= 17730918] = max(all.hor.suplications$duplication_id) + 1
  
  
  
}
### ### ### ### ### 




duplication.summary = data.frame(duplication.id = vector(mode = "numeric", length = max(all.hor.suplications$duplication_id)),
                                 accession = vector(mode = "numeric", length = max(all.hor.suplications$duplication_id)),
                                 chromosome = vector(mode = "numeric", length = max(all.hor.suplications$duplication_id)),
                                 min.startA = vector(mode = "numeric", length = max(all.hor.suplications$duplication_id)),
                                 min.startB = vector(mode = "numeric", length = max(all.hor.suplications$duplication_id)),
                                 max.endA = vector(mode = "numeric", length = max(all.hor.suplications$duplication_id)),
                                 max.endB = vector(mode = "numeric", length = max(all.hor.suplications$duplication_id)),
                                 direction = vector(mode = "numeric", length = max(all.hor.suplications$duplication_id)),
                                 block.size.total = vector(mode = "numeric", length = max(all.hor.suplications$duplication_id)),
                                 estimated.size = vector(mode = "numeric", length = max(all.hor.suplications$duplication_id)),
                                 no.of.HORs = vector(mode = "numeric", length = max(all.hor.suplications$duplication_id)))
for(i in 1 : nrow(duplication.summary))
{
  hor.temp = all.hor.suplications[all.hor.suplications$duplication_id == i,]
  
  duplication.summary$duplication.id[i] = i
  duplication.summary$accession[i] = hor.temp$accession[1]
  duplication.summary$chromosome[i] = strsplit(hor.temp$chrA[1], split = "")[[1]][length(strsplit(hor.temp$chrA[1], split = "")[[1]])]
  duplication.summary$min.startA[i] = min(hor.temp$start.A.bp)
  duplication.summary$min.startB[i] = min(hor.temp$start.B.bp)
  duplication.summary$max.endA[i] = max(hor.temp$end.A.bp)
  duplication.summary$max.endB[i] = max(hor.temp$end.B.bp)
  duplication.summary$direction[i] = hor.temp$direction.1.para_2.perp.[1]
  duplication.summary$block.size.total[i] = sum(hor.temp$block.size.in.units)
  duplication.summary$no.of.HORs[i] = nrow(hor.temp)
  duplication.summary$estimated.size[i] = ave(duplication.summary$max.endA[i] - duplication.summary$min.startA[i], duplication.summary$max.endB[i] - duplication.summary$min.startB[i])
  
}


duplication.summary = duplication.summary[duplication.summary$no.of.HORs >= 50,]


write.csv(x = duplication.summary, file = "duplication.summary.csv")
duplication.summary = read.csv(file = "duplication.summary.csv")


TE = read.table("66Atha_fulllength+soloLTRs.matrix+phylo+hmm+cengap", header = TRUE)
HOR.files.accession %in% TE$accession_fasta




#plot again
{
  pdf(file = "t0c5 HOR 0 variant extra lines only containing.pdf", onefile = TRUE, width = 14)
  par(mfrow = c(1,2))
  for(i in 1 : length(HOR.files))
  {
    chromosome = HOR.files.chr[i]
    assembly = HOR.files.accession[i]
    print(paste(i, ": working on chromosome ", chromosome, " from ", assembly, sep = ""))
    hors = read.csv(HOR.files[i])
    #repeats = read.csv(file = repeat.files[which(repeat.files.accession == assembly & repeat.files.chr == chromosome)])
    
    hors = hors[!(hors$start.A.bp  %in%  boxplot.stats(hors$start.A.bp)$out),]
    hors = hors[!(hors$start.B.bp  %in%  boxplot.stats(hors$start.B.bp)$out),]
    
    hors = hors[hors$total_variant <= 0,]
    hors = hors[hors$block.size.in.units >= 5,]
    hors.save = hors
    
    if(assembly %in% duplication.summary$accession & chromosome %in% duplication.summary$chromosome[duplication.summary$accession == assembly])
    {
      
      if(nrow(hors) > 0)
      {
        plot(NA, main = paste("chr ", chromosome, " from ", assembly, sep = ""),
             xlim = c(min(min(hors$start.A.bp), min(hors$start.B.bp)), max(max(hors$end.A.bp), max(hors$end.B.bp))),
             ylim = c(min(min(hors$start.A.bp), min(hors$start.B.bp)), max(max(hors$end.A.bp), max(hors$end.B.bp))))
        for(j in 1 : nrow(hors))
        {
          lines(c(hors$start.A.bp[j], hors$end.A.bp[j]), c(hors$start.B.bp[j], hors$end.B.bp[j]), lwd = 3)
        }
      }
      
      temp.duplications = duplication.summary[(duplication.summary$accession == assembly & duplication.summary$chromosome == chromosome),]
      if(nrow(temp.duplications) > 0)
      {
        for(j in 1 : nrow(temp.duplications))
        {
          if(temp.duplications$direction[j] == 2)
          {
            lines(c(temp.duplications$min.startA[j], temp.duplications$max.endA[j]), c(temp.duplications$max.endB[j], temp.duplications$min.startB[j]), lwd = 1.5, col = "red")
            
          } else
          {
            lines(c(temp.duplications$min.startA[j], temp.duplications$max.endA[j]), c(temp.duplications$min.startB[j], temp.duplications$max.endB[j]), lwd = 1.5, col = "red")
            
          }
        }
      }
      
      TE.temp = TE[TE$accession_fasta == assembly & TE$chromosome == paste("Chr", chromosome, sep = ""),]
      TE.temp = TE.temp[TE.temp$genome_left_coord_FL >= min(min(hors$start.A.bp), min(hors$start.B.bp)),]
      TE.temp = TE.temp[TE.temp$genome_right_coord_FL <= max(max(hors$end.A.bp), max(hors$end.B.bp)),]
      if(nrow(TE.temp) > 0)
      {
        for(j in 1 : nrow(TE.temp))
        {
          if(TE.temp$quality[j] == "intact")
          {
            abline(h = TE.temp$genome_left_coord_FL[j], v = TE.temp$genome_left_coord_FL[j], col = "blue", lty = 3)
          } else
          {
            abline(h = TE.temp$genome_left_coord_FL[j], v = TE.temp$genome_left_coord_FL[j], col = "red", lty = 3)
          }
        }
      }
      
      hist(abs(hors.save$start_A - hors.save$start_B), breaks = 50, main = "distances between blocks")
      
    }
    
    
    
  }
  dev.off()
}


duplication.summary$chromosome = as.numeric(duplication.summary$chromosome)
duplication.summary$est.identity = (duplication.summary$block.size.total * 178) / duplication.summary$estimated.size


pdf(file = "duplications 66 chromosome vs size.pdf")
stripchart((duplication.summary$estimated.size/1000) ~ duplication.summary$chromosome, pch = 16, method = "jitter",
           xlab = "chromosome number",
           ylab = "duplication size, kb", vertical=TRUE,
           main = "duplications: chromosome vs size")
dev.off()

plot(duplication.summary$est.identity, duplication.summary$estimated.size)





### ### ### ### ### ### 
### plot again with the other half of the plot having t5c2 hors
### ### ### ### ### ### 



hor.more.files = list.files(path = "C:/Users/wlodz/Desktop/panC_analysis_athome/detect_large_duplications/t5c3", pattern = "HORs_", full.names = TRUE)

#plot again
{
  #pdf(file = "t0c5 HOR and t5c3 extra lines all.pdf", onefile = TRUE, width = 14)

  for(i in 1 : length(HOR.files))
  {
    #i = 314
    pdf(file = paste("plots/", assembly, "_", chromosome, "_t0c5 HOR and t5c3 extra lines all.pdf", sep = ""))
    
    chromosome = HOR.files.chr[i]
    assembly = HOR.files.accession[i]
    #png(filename = paste("plots/", assembly, "_", chromosome, "_t0c5 HOR and t5c3 extra lines all.png", sep = ""), width = 2000, height = 2000)
    print(paste(i, ": working on chromosome ", chromosome, " from ", assembly, sep = ""))
    hors = read.csv(HOR.files[i])
    hors.more = read.csv(hor.more.files[i])
    #repeats = read.csv(file = repeat.files[which(repeat.files.accession == assembly & repeat.files.chr == chromosome)])
    
    #hors = hors[!(hors$start.A.bp  %in%  boxplot.stats(hors$start.A.bp)$out),]
    #hors = hors[!(hors$start.B.bp  %in%  boxplot.stats(hors$start.B.bp)$out),]
    #hors.more = hors.more[!(hors.more$start.A.bp  %in%  boxplot.stats(hors.more$start.A.bp)$out),]
    #hors.more = hors.more[!(hors.more$start.B.bp  %in%  boxplot.stats(hors.more$start.B.bp)$out),]
    
    hors.more = hors.more[sample(nrow(hors.more), size = round(nrow(hors.more))/10),]
    
    hors = hors[hors$total_variant <= 0,]
    hors = hors[hors$block.size.in.units >= 5,]
    
    if(nrow(hors) > 0 & nrow(hors.more) > 0)
    {
      limsBoth = c(c(min(min(hors.more$start.A.bp), min(hors.more$start.B.bp)), max(max(hors.more$end.A.bp), max(hors.more$end.B.bp))))
      #limsBoth = c(12328778, 18525923)
      
      
      plot(NA, main = paste("chr ", chromosome, " from ", assembly, sep = ""),
           xlim = limsBoth,
           ylim = limsBoth)
      #points(hors.more$start.B.bp, hors.more$start.A.bp, pch = 20, col = "grey", cex = 0.5)
      for(j in 1 : nrow(hors.more))
      {
        lines(c(hors.more$start.B.bp[j], hors.more$end.B.bp[j]), c(hors.more$start.A.bp[j], hors.more$end.A.bp[j]), lwd = 3, col = "grey")
      }
      for(j in 1 : nrow(hors))
      {
        lines(c(hors$start.A.bp[j], hors$end.A.bp[j]), c(hors$start.B.bp[j], hors$end.B.bp[j]), lwd = 3)
      }
      
      
    }
    temp.duplications = duplication.summary[(duplication.summary$accession == assembly & duplication.summary$chromosome == chromosome),]
    if(nrow(temp.duplications) > 0)
    {
      for(j in 1 : nrow(temp.duplications))
      {
        if(temp.duplications$direction[j] == 2)
        {
          lines(c(temp.duplications$min.startA[j], temp.duplications$max.endA[j]), c(temp.duplications$max.endB[j], temp.duplications$min.startB[j]), lwd = 1.5, col = "red")
          
        } else
        {
          lines(c(temp.duplications$min.startA[j], temp.duplications$max.endA[j]), c(temp.duplications$min.startB[j], temp.duplications$max.endB[j]), lwd = 1.5, col = "red")
          
        }
      }
    }
    TE.temp = TE[TE$accession_fasta == assembly & TE$chromosome == paste("Chr", chromosome, sep = ""),]
    TE.temp = TE.temp[TE.temp$genome_left_coord_FL >= min(min(hors.more$start.A.bp), min(hors.more$start.B.bp)),]
    TE.temp = TE.temp[TE.temp$genome_right_coord_FL <= max(max(hors.more$end.A.bp), max(hors.more$end.B.bp)),]
    if(nrow(TE.temp) > 0)
    {
      for(j in 1 : nrow(TE.temp))
      {
        if(TE.temp$quality[j] == "intact")
        {
          abline(h = TE.temp$genome_left_coord_FL[j], v = TE.temp$genome_left_coord_FL[j], col = "blue", lty = 3)
        } else
        {
          abline(h = TE.temp$genome_left_coord_FL[j], v = TE.temp$genome_left_coord_FL[j], col = "red", lty = 3)
        }
      }
    }
    dev.off()
  }
  #dev.off()
}


#######################################
################# OLD #################
#######################################

#check for block distance normality in windows of f repeats

repeats_f = 500


for(i in 1 : length(HOR.files))
{
  repeats = read.csv(file = repeat.files[i])
  hors = read.csv(file = HOR.files[i])
  hors = hors[order(hors$start_A),]
  hors$block_distance = hors$start.B.bp - hors$start.A.bp
  quantiles.variant.score = quantile(hors$total_variant, seq(0,1,0.01))
  hors$in.quantile = FALSE
  hors$in.quantile[hors$total_variant < quantiles.variant.score[2]] = TRUE
  
  
  
  block_distance_windows_normality = vector(mode = "numeric", length = ceiling(nrow(repeats)/repeats_f))
  
  limX = c(min(hors$start.A.bp), max(hors$start.A.bp))
  limY = c(min(hors$start.B.bp), max(hors$start.B.bp))
  
  png(filename = paste("plots/histograms_", repeats$assembly[1], "_", repeats$chromosome[i], "_", repeats_f, ".png"), 
      width = 2000, height = 2000)
  plot(hors$start.A.bp, hors$start.B.bp, type = "p", pch = 20, col = "black", xlim = limX, ylim = limY)
  par(new = TRUE)
  plot(hors$start.A.bp[hors$in.quantile], hors$start.B.bp[hors$in.quantile], type = "p", pch = 20, col = "blue", xlim = limX, ylim = limY, xlab = "", ylab = "")
  dev.off()
  
  pdf(file = paste("plots/histograms_", repeats$assembly[1], "_", repeats$chromosome[i], "_", repeats_f, ".pdf"), onefile = TRUE)
  plot(hors$start.B.bp[hors$in.quantile] - hors$start.A.bp[hors$in.quantile], type = "p", pch = 20, col = "blue")
  hist(hors$start.B.bp[hors$in.quantile] - hors$start.A.bp[hors$in.quantile], breaks = 100)
  
  for(j in 1 : length(block_distance_windows_normality))
  {
    print(paste(i, ": ",j, "/", length(block_distance_windows_normality)))
    
    
    hor_subset = hors[((hors$start_A > (j-1)*repeats_f) & (hors$start_A <= (j)*repeats_f)) | ((hors$start_B > (j-1)*repeats_f) & (hors$start_B <= (j)*repeats_f)),]
    hor_subset = hor_subset[hor_subset$in.quantile,]
      
      if(nrow(hor_subset) >= 5000)
    {
      hor_subset = hor_subset[sample(nrow(hor_subset), 5000),]
    }
    if(nrow(hor_subset) >= 3)
    {
      #block_distance_windows_normality[j] = shapiro.test(hor_subset$block_distance)$p.value #this is almost never normal
      hist(hor_subset$block_distance, breaks = 100, xlim = c(0, max(hors$block_distance)), ylim = c(0,5))
    }
  }
  dev.off()
}

hors.all = read.csv(file = HOR.files[1])
hors.all = hors[0,]
for(i in 1 : length(HOR.files))
{
  print(i)
  hors = read.csv(file = HOR.files[i])
  
  quantiles.variant.score = quantile(hors$total_variant, seq(0,1,0.01))
  hors$in.quantile = FALSE
  hors$in.quantile[hors$total_variant < quantiles.variant.score[2]] = TRUE
  
  hors.all = rbind(hors.all, hors[hors$in.quantile,])
}

hors.all$block_distance = hors.all$start.B.bp - hors.all$start.A.bp


pdf(file = "all.young.hors.distance.histogram.pdf", onefile = TRUE)
hist(hors.all$block_distance, breaks = 100)
hist(hors.all$block.size.in.units, breaks = 100)
dev.off()


library(tidyverse)

a=hist(hors.all$block_distance, breaks = 1000)

df <- tibble(t = 1:length(a$counts), y = a$counts, sensor = 'sensor1')

fit <- nls(y ~ SSasymp(t, yf, y0, log_alpha), data = df)







b = data.frame(t = a$breaks[1: (length(a$breaks) - 1)], y = a$counts)
b$data = "data1"


c <- tibble(t = a$breaks[1: (length(a$breaks) - 1)], y = a$counts, sensor = 'sensor1')

fit <- nls(y ~ SSasymp(t, yf, y0, log_alpha), data = sensor1)

















