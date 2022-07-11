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

setwd("C:/Users/wlodz/Desktop/panC_analysis_athome/athila_insertion")

skip.analysis = TRUE

if(!skip.analysis)
{
  
  assemblies = list.files(path = "C:/Users/wlodz/Desktop/panC_analysis_athome/assemblies", pattern = ".f", full.names = T)
  assemblies.names = list.files(path = "C:/Users/wlodz/Desktop/panC_analysis_athome/assemblies", pattern = ".f", full.names = F)
  
  
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
  
  #TEo = read.table("C:/Users/wlodz/Desktop/panC_analysis_athome/detect_large_duplications/66Atha_fulllength+soloLTRs.matrix+phylo+hmm+cengap", header = TRUE)
  TE = read.table("66Atha_fulllength+soloLTRs.matrix+phylo+hmm+cengap+solofam_new", header = TRUE)
  TE$up20 = ""
  TE$down20 = ""
  
  TE$cen180up = ""
  TE$cen180down = ""
  
  TE$strand.up = ""
  TE$strand.down = ""
  
  
  #repeat.files.accession %in% TE$accession_fasta
  #assemblies.names %in% TE$accession_fasta
  
  
  #for each athila
  
  for(i in 1 : length(assemblies))
  {
    assembly.name = assemblies.names[i]
    assembly = read.fasta(file = assemblies[i])
    for(j in 1 : 5)
    {
      chromosome = j
      repeats = read.csv(repeat.files[repeat.files.accession == assembly.name & repeat.files.chr == as.character(chromosome)])
      
      athilas.todo = which(TE$accession_fasta == assembly.name & TE$chromosome == paste("Chr", chromosome, sep = ""))
      
      for(k in athilas.todo)
      {
        print(paste(i, j, ": Athila ", k, "/", nrow(TE), sep = ""))
        
        start.c = TE$genome_left_coord_FL[k]
        end.c = TE$genome_right_coord_FL[k]
        
        if(start.c <= 0)
        {
          start.c = 1
        }
        if(end.c > length(assembly[[chromosome]]))
        {
          end.c = length(assembly[[chromosome]])
        }
        
        TE$up20[k] = paste(assembly[[chromosome]][(start.c - 20) : (start.c - 1)], collapse = "")
        TE$down20[k] = paste(assembly[[chromosome]][(end.c + 1) : (end.c + 20)], collapse = "")
        
        if(length(which(repeats$start > end.c)) > 0)
        {
          first.cen180.up = min(which(repeats$start > end.c))
          TE$cen180up[k] = repeats$sequence_strand_adjusted[first.cen180.up]
          TE$strand.up[k] = repeats$strand[first.cen180.up]
        } else
        {
          TE$cen180up[k] = ""
          TE$strand.up[k] = ""
        }
        
        if(length(which(repeats$start < end.c)) > 0)
        {
          last.cen180.down = max(which(repeats$start < end.c))
          TE$cen180down[k] = repeats$sequence_strand_adjusted[first.cen180.up]
          TE$strand.down[k] = repeats$strand[last.cen180.down]
        } else
        {
          TE$cen180down[k] = ""
          TE$strand.down[k] = ""
        }
      }
    }
  }
  
  write.csv(x = TE, file = "TE_with20bp.csv")
  TE = read.csv(file = "TE_with20bp.csv", header = T)
  
  #see which athilas have repeats on opposing strands around, meaning there is a switch in directions
  {
    
    TE$strand.up[TE$strand.up == "+"] = 1
    TE$strand.up[TE$strand.up == "-"] = 2
    TE$strand.down[TE$strand.down == "+"] = 1
    TE$strand.down[TE$strand.down == "-"] = 2
    
    TE$strand.down = as.numeric(TE$strand.down)
    TE$strand.up = as.numeric(TE$strand.up)
    
    TE$is.same.strand.cen180[TE$strand.down == 1 & TE$strand.up == 2] = F
    TE$is.same.strand.cen180[TE$strand.down == 2 & TE$strand.up == 1] = F
    TE$is.same.strand.cen180[TE$strand.down == 1 & TE$strand.up == 1] = T
    TE$is.same.strand.cen180[TE$strand.down == 2 & TE$strand.up == 2] = T
    
    length(which(TE$is.same.strand.cen180 == T))
    #8608
    length(which(TE$is.same.strand.cen180 == F))
    #6418
    
    TE$strand.up[TE$strand.up == 1] = "+"
    TE$strand.up[TE$strand.up == 2] = "-"
    TE$strand.down[TE$strand.down == 1] = "+"
    TE$strand.down[TE$strand.down == 2] = "-"
  }
  
  #align and score the alignments between 20bp flanks and next closest repeat
  
  TE$position.on.180.up = 0
  TE$position.on.180.down = 0
  TE$score.up = NA
  TE$score.down = NA
  
  for(i in 1 : nrow(TE))
  {
    print(i)
    
    if(!is.na(TE$strand.up[i]) & !is.na(TE$cen180up[i]))
    {
      seq20 = toupper(TE$up20[i])
      cen180 = TE$cen180up[i]
      strand = TE$strand.up[i]
      if(strand == "-")
      {
        alignment = pairwiseAlignment(pattern = revCompString(seq20), cen180, type = "global-local", gapOpening = 10)
        TE$score.up[i] = score(alignment)
        TE$position.on.180.up[i] = start(subject(alignment))
      } else
      {
        alignment = pairwiseAlignment(pattern = seq20, cen180, type = "global-local", gapOpening = 10)
        TE$score.up[i] = score(alignment)
        TE$position.on.180.up[i] = start(subject(alignment))
      }
    }
    
    if(!is.na(TE$strand.down[i]) & !is.na(TE$cen180down[i]))
    {
      seq20 = toupper(TE$down20[i])
      cen180 = TE$cen180down[i]
      strand = TE$strand.down[i]
      if(strand == "-")
      {
        alignment = pairwiseAlignment(pattern = revCompString(seq20), cen180, type = "global-local", gapOpening = 10)
        TE$score.down[i] = score(alignment)
        TE$position.on.180.down[i] = start(subject(alignment))
      } else
      {
        alignment = pairwiseAlignment(pattern = seq20, cen180, type = "global-local", gapOpening = 10)
        TE$score.down[i] = score(alignment)
        TE$position.on.180.down[i] = start(subject(alignment))
      }
    }
  }
  
  #align and score the alignments between 20bp flanks and cen180 consensus (expanded by 20 bp both sides)
  
  TE$position.on.180.up.consensus = 0
  TE$position.on.180.down.consensus = 0
  TE$score.up.consensus = NA
  TE$score.down.consensus = NA
  cen180consensus = "AGTATAAGAACTTAAACCGCAACCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAAGACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAGAAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTGAGTATAAGAACTTAAACCGC"
  cen180 = cen180consensus
  
  for(i in 1 : nrow(TE))
  {
    if(TE$centromere[i] == "in")
    {
      print(i)
      
      if(!is.na(TE$strand.up[i]) & !is.na(TE$cen180up[i]))
      {
        seq20 = toupper(TE$up20[i])
        #cen180 = TE$cen180up[i]
        strand = TE$strand.up[i]
        if(strand == "-")
        {
          alignment = pairwiseAlignment(pattern = revCompString(seq20), cen180, type = "global-local", gapOpening = 10)
          TE$score.up.consensus[i] = score(alignment)
          TE$position.on.180.up.consensus[i] = start(subject(alignment))
        } else
        {
          alignment = pairwiseAlignment(pattern = seq20, cen180, type = "global-local", gapOpening = 10)
          TE$score.up.consensus[i] = score(alignment)
          TE$position.on.180.up.consensus[i] = start(subject(alignment))
        }
      }
      
      if(!is.na(TE$strand.down[i]) & !is.na(TE$cen180down[i]))
      {
        seq20 = toupper(TE$down20[i])
        #cen180 = TE$cen180down[i]
        strand = TE$strand.down[i]
        if(strand == "-")
        {
          alignment = pairwiseAlignment(pattern = revCompString(seq20), cen180, type = "global-local", gapOpening = 10)
          TE$score.down.consensus[i] = score(alignment)
          TE$position.on.180.down.consensus[i] = start(subject(alignment))
        } else
        {
          alignment = pairwiseAlignment(pattern = seq20, cen180, type = "global-local", gapOpening = 10)
          TE$score.down.consensus[i] = score(alignment)
          TE$position.on.180.down.consensus[i] = start(subject(alignment))
        }
      }
    }
  }
  TE$position.on.180.up.consensus[TE$position.on.180.up.consensus > 177] = TE$position.on.180.up.consensus[TE$position.on.180.up.consensus > 177] - 177
  TE$position.on.180.down.consensus[TE$position.on.180.down.consensus > 177] = TE$position.on.180.down.consensus[TE$position.on.180.down.consensus > 177] - 177
  
  
  write.csv(x = TE, file = "TE_with_flanks_scored.csv")
}


TE = read.csv(file = "TE_with_flanks_scored.csv")

TE = TE[TE$centromere == "in",]
TE = columnRemove(TE, c("X"))

hist(c(TE$score.up.consensus,TE$score.down.consensus), breaks = 80)


### ### ###
#uncomment the 2 below to use alignments with closest cen180 instead of consensus cen180, it makes no difference though
### ### ###

#TE$position.on.180.up.consensus = TE$position.on.180.up
#TE$position.on.180.down.consensus = TE$position.on.180.down


#correct the position by adding or subtracting bp based on the strand for both up and down
TE$strand.down[is.na(TE$strand.down)] = "o" #easier to handle than NAs
TE$strand.up[is.na(TE$strand.up)] = "o" 

TE$position.on.180.up.consensus[TE$position.on.180.up.consensus == 0] = NA     #position cannot be 0 as it could be a valid position, need to keep as NA
TE$position.on.180.down.consensus[TE$position.on.180.down.consensus == 0] = NA

TE$position.on.180.up.consensus.adjusted[TE$strand.up == "+"] = TE$position.on.180.up.consensus[TE$strand.up == "+"] + 19 #end of the alignment counts as the break site, so gotta add 19
TE$position.on.180.up.consensus.adjusted[TE$strand.up == "-"] = TE$position.on.180.up.consensus[TE$strand.up == "-"]
TE$position.on.180.up.consensus.adjusted[TE$position.on.180.up.consensus.adjusted > 177 & TE$strand.up == "+"] = TE$position.on.180.up.consensus.adjusted[TE$position.on.180.up.consensus.adjusted > 177 & TE$strand.up == "+"] - 177 #some are getting longer than consensus size, need to be wrapped back

TE$position.on.180.down.consensus.adjusted[TE$strand.down == "+" & !is.na(TE$position.on.180.down.consensus)] = TE$position.on.180.down.consensus[TE$strand.down == "+" & !is.na(TE$position.on.180.down.consensus)] 
TE$position.on.180.down.consensus.adjusted[TE$strand.down == "-" & !is.na(TE$position.on.180.down.consensus)] = TE$position.on.180.down.consensus[TE$strand.down == "-" & !is.na(TE$position.on.180.down.consensus)] + 19
TE$position.on.180.down.consensus.adjusted[TE$position.on.180.down.consensus.adjusted > 177 & TE$strand.down == "-" & !is.na(TE$position.on.180.down.consensus)] = TE$position.on.180.down.consensus.adjusted[TE$position.on.180.down.consensus.adjusted > 177 & TE$strand.down == "-" & !is.na(TE$position.on.180.down.consensus)] - 177


#calculate distance between them

which.are.same.strand.plus = (!is.na(TE$position.on.180.up.consensus.adjusted) & 
                                !is.na(TE$position.on.180.down.consensus.adjusted) & 
                                TE$strand.up == TE$strand.down & 
                                TE$strand.up == "+")
TE$distance.to.cut.sum[which.are.same.strand.plus] = TE$position.on.180.down.consensus.adjusted[which.are.same.strand.plus] - TE$position.on.180.up.consensus.adjusted[which.are.same.strand.plus]


which.are.same.strand.minus = (!is.na(TE$position.on.180.up.consensus.adjusted) & 
                                 !is.na(TE$position.on.180.down.consensus.adjusted) & 
                                 TE$strand.up == TE$strand.down & 
                                 TE$strand.up == "-")
TE$distance.to.cut.sum[which.are.same.strand.minus] = TE$position.on.180.up.consensus.adjusted[which.are.same.strand.minus] - TE$position.on.180.down.consensus.adjusted[which.are.same.strand.minus]

save.again = FALSE
if(save.again)
{
  write.csv(x = TE, file = "TE_with_flanks_scored_adjusted.csv")
}


### ### ###
# plot
### ### ###

#set1: all that pair

TEtemp = TE[!is.na(TE$position.on.180.down.consensus.adjusted) & !is.na(TE$position.on.180.up.consensus.adjusted),]

hist(TEtemp$distance.to.cut.sum[!is.na(TEtemp$distance.to.cut.sum)], breaks = 400, xlim = c(-200,200))

#set2: all that pair and high score

TEtemp = TE[!is.na(TE$position.on.180.down.consensus.adjusted) & 
              !is.na(TE$position.on.180.up.consensus.adjusted) &
              (TE$score.up.consensus >= 0) & !is.na(TE$score.up.consensus) &
              (TE$score.down.consensus >= 0) & !is.na(TE$score.down.consensus) &
              TE$centromere == "in",]

hist(TEtemp$distance.to.cut.sum[!is.na(TEtemp$distance.to.cut.sum)], breaks = 400, xlim = c(-200,200))

pdf(file = "histogram.insertion.site.size.difference.pdf")
hist(TEtemp$distance.to.cut.sum[!is.na(TEtemp$distance.to.cut.sum)] - 1, breaks = 400, xlim = c(-200,200),
     main = "histogram of distances between flanking regions")

hist(c(TEtemp$position.on.180.down.consensus.adjusted[!is.na(TEtemp$position.on.180.down.consensus.adjusted)], 
       TEtemp$position.on.180.up.consensus.adjusted[!is.na(TEtemp$position.on.180.up.consensus.adjusted)]), 
     breaks = 300, 
     xlim = c(0,178),
     main = "histogram of insertion sites along cen180 consensus")
dev.off()

#set3: pair, high score, same strand + or -

TEtemp = TE[!is.na(TE$position.on.180.down.consensus.adjusted) & 
              !is.na(TE$position.on.180.up.consensus.adjusted) &
              (TE$score.up.consensus >= 0) & !is.na(TE$score.up.consensus) &
              (TE$score.down.consensus >= 0) & !is.na(TE$score.down.consensus) &
              TE$strand.up == "+",]
hist(TEtemp$distance.to.cut.sum[!is.na(TEtemp$distance.to.cut.sum)], breaks = 400, xlim = c(-200,200))


TEtemp = TE[!is.na(TE$position.on.180.down.consensus.adjusted) & 
              !is.na(TE$position.on.180.up.consensus.adjusted) &
              (TE$score.up.consensus >= 0) & !is.na(TE$score.up.consensus) &
              (TE$score.down.consensus >= 0) & !is.na(TE$score.down.consensus) &
              TE$strand.up == "-",]
hist(TEtemp$distance.to.cut.sum[!is.na(TEtemp$distance.to.cut.sum)], breaks = 400, xlim = c(-200,200))


#set4: different TE chromosomes
pdf(file = paste("Chromosomes separately histogram.insertion.site.size.difference.pdf", sep = ""))
par(mfrow = c(2,1))
for(i in 1 : 5)
{
  print(i)
  TEtemp = TE[!is.na(TE$position.on.180.down.consensus.adjusted) & 
                !is.na(TE$position.on.180.up.consensus.adjusted) &
                (TE$score.up.consensus >= 0) & !is.na(TE$score.up.consensus) &
                (TE$score.down.consensus >= 0) & !is.na(TE$score.down.consensus) &
                TE$chromosome == paste("Chr", i, sep = ""),]
  
  if(nrow(TEtemp) > 0)
  {
    
    hist(TEtemp$distance.to.cut.sum[!is.na(TEtemp$distance.to.cut.sum)] - 1, breaks = 400, xlim = c(-200,200),
         main =  paste(paste("Chr", i, sep = "")," histogram of distances between flanking regions", sep = ""))
    
    hist(c(TEtemp$position.on.180.down.consensus.adjusted[!is.na(TEtemp$position.on.180.down.consensus.adjusted)], 
           TEtemp$position.on.180.up.consensus.adjusted[!is.na(TEtemp$position.on.180.up.consensus.adjusted)]), 
         breaks = 300, 
         xlim = c(0,178),
         main =  paste(paste("Chr", i, sep = "")," histogram of insertion sites along cen180 consensus", sep = ""))
    
  } else
  {
    print("did not find enough TEs with specified conditions")
  }
}
dev.off()

#set5: different TE families in centromere

pdf(file = paste("families.in.centromere.histograms.insertion.sites.pdf", sep = ""))
par(mfrow = c(2,1))
for(i in 1 : length(unique(TE$combo_fam)))
{
  print(i)
  TEtemp = TE[!is.na(TE$position.on.180.down.consensus.adjusted) & 
                !is.na(TE$position.on.180.up.consensus.adjusted) &
                (TE$score.up.consensus >= 0) & !is.na(TE$score.up.consensus) &
                (TE$score.down.consensus >= 0) & !is.na(TE$score.down.consensus) &
                TE$combo_fam == unique(TE$combo_fam)[i] &
                TE$centromere == "in",]
  
  if(nrow(TEtemp) > 0)
  {
    hist(TEtemp$distance.to.cut.sum[!is.na(TEtemp$distance.to.cut.sum)] - 1, breaks = 400, xlim = c(-200,200),
         main =  paste("family ", unique(TE$combo_fam)[i]," histogram of distances between flanking regions", sep = ""))
    
    hist(c(TEtemp$position.on.180.down.consensus.adjusted[!is.na(TEtemp$position.on.180.down.consensus.adjusted)], 
           TEtemp$position.on.180.up.consensus.adjusted[!is.na(TEtemp$position.on.180.up.consensus.adjusted)]), 
         breaks = 300, 
         xlim = c(0,178),
         main =  paste("family ", unique(TE$combo_fam)[i]," histogram of insertion sites along cen180 consensus", sep = ""))
    
  } else
  {
    print("did not find enough TEs with specified conditions")
  }
  
}
dev.off()


pdf(file = paste("intact.families.in.centromere.histograms.insertion.sites.pdf", sep = ""))
par(mfrow = c(2,1))
for(i in 1 : length(unique(TE$combo_fam)))
{
  print(i)
  TEtemp = TE[!is.na(TE$position.on.180.down.consensus.adjusted) & 
                !is.na(TE$position.on.180.up.consensus.adjusted) &
                (TE$score.up.consensus >= 0) & !is.na(TE$score.up.consensus) &
                (TE$score.down.consensus >= 0) & !is.na(TE$score.down.consensus) &
                TE$combo_fam == unique(TE$combo_fam)[i] &
                TE$centromere == "in" &
                TE$quality == "intact",]
  
  if(nrow(TEtemp) > 0)
  {
    hist(TEtemp$distance.to.cut.sum[!is.na(TEtemp$distance.to.cut.sum)] - 1, breaks = 400, xlim = c(-200,200),
         main =  paste("family ", unique(TE$combo_fam)[i]," histogram of distances between flanking regions", sep = ""))
    
    hist(c(TEtemp$position.on.180.down.consensus.adjusted[!is.na(TEtemp$position.on.180.down.consensus.adjusted)], 
           TEtemp$position.on.180.up.consensus.adjusted[!is.na(TEtemp$position.on.180.up.consensus.adjusted)]), 
         breaks = 300, 
         xlim = c(0,178),
         main =  paste("family ", unique(TE$combo_fam)[i]," histogram of insertion sites along cen180 consensus", sep = ""))
    
  } else
  {
    print("did not find enough TEs with specified conditions")
  }
  
}
dev.off()

pdf(file = paste("solo.families.in.centromere.histograms.insertion.sites.pdf", sep = ""))
par(mfrow = c(2,1))
for(i in 1 : length(unique(TE$combo_fam)))
{
  print(i)
  TEtemp = TE[!is.na(TE$position.on.180.down.consensus.adjusted) & 
                !is.na(TE$position.on.180.up.consensus.adjusted) &
                (TE$score.up.consensus >= 0) & !is.na(TE$score.up.consensus) &
                (TE$score.down.consensus >= 0) & !is.na(TE$score.down.consensus) &
                TE$combo_fam == unique(TE$combo_fam)[i] &
                TE$centromere == "in" &
                TE$quality == "solo",]
  
  if(nrow(TEtemp) > 0)
  {
    hist(TEtemp$distance.to.cut.sum[!is.na(TEtemp$distance.to.cut.sum)] - 1, breaks = 400, xlim = c(-200,200),
         main =  paste("family ", unique(TE$combo_fam)[i]," histogram of distances between flanking regions", sep = ""))
    
    hist(c(TEtemp$position.on.180.down.consensus.adjusted[!is.na(TEtemp$position.on.180.down.consensus.adjusted)], 
           TEtemp$position.on.180.up.consensus.adjusted[!is.na(TEtemp$position.on.180.up.consensus.adjusted)]), 
         breaks = 300, 
         xlim = c(0,178),
         main =  paste("family ", unique(TE$combo_fam)[i]," histogram of insertion sites along cen180 consensus", sep = ""))
    
  } else
  {
    print("did not find enough TEs with specified conditions")
  }
  
}
dev.off()




