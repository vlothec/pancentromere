library(stringr)
library(base)
library(msa)
library(Biostrings)
library(seqinr)


setwd("C:/Users/wlodz/Desktop/panC_analysis_athome/non_178_location")



repeats = read.csv(file = "C:/Users/wlodz/Desktop/panC_analysis_athome/normality/all_repeats_66_consensus_HOR.csv")

relative.positions = NULL

for(i in 1 : length(unique(repeats$assembly)))
{
  repeat.sample = repeats[repeats$assembly == unique(repeats$assembly)[i],]
  for(j in 1)# : 5)
  {
    repeat.sample.chr = repeat.sample[repeat.sample$chromosome == unique(repeat.sample$chromosome)[j],]
    
    relative.positions = c(relative.positions, which(repeat.sample.chr$len != 178)/nrow(repeat.sample.chr))
    
    print(paste(((i-1)*5) + j, sep = ""))
  }
}

pdf(file = "chr1_relative_pos_non_178_repeats.pdf", width = 10, height = 5)
hist(relative.positions)
dev.off()



relative.positions = NULL

for(i in 1 : length(unique(repeats$assembly)))
{
  repeat.sample = repeats[repeats$assembly == unique(repeats$assembly)[i],]
  for(j in 1)# : 5)
  {
    repeat.sample.chr = repeat.sample[repeat.sample$chromosome == unique(repeat.sample$chromosome)[j],]
    
    relative.positions = c(relative.positions, which(repeat.sample.chr$HORlengthsSum > 500)/nrow(repeat.sample.chr))
    
    print(paste(((i-1)*5) + j, sep = ""))
  }
}

pdf(file = "chr1_relative_pos_HOR_500plus_repeats.pdf", width = 10, height = 5)
hist(relative.positions)
dev.off()



pdf(file = "chr1_repeats HOR values 178 vs non 178.pdf", onefile = TRUE)
hist(repeats$HORlengthsSum[repeats$len != 178], breaks = 20000, xlim = c(50,3000), ylim = c(0,30000))
hist(repeats$HORlengthsSum[repeats$len == 178], breaks = 20000, xlim = c(05,3000), ylim = c(0,30000))
dev.off()


repeats$len = as.factor(repeats$len)

repeats$HORlengthsSum = repeats$HORlengthsSum + 1
p = ggplot(repeats[repeats$length > 160 & repeats$length < 200,], aes(x = len, y = log2(HORlengthsSum))) + 
  geom_violin() + 
  geom_boxplot(width=0.1) +
  geom_hline(yintercept = mean(log2(repeats$HORlengthsSum)))
ggsave(filename = "chr1_log2HOR violin sizes more.pdf", width = 30, height = 6)
repeats$HORlengthsSum = repeats$HORlengthsSum - 1


repeats$edit.distance = repeats$edit.distance + 1
p = ggplot(repeats[repeats$length > 160 & repeats$length < 200,], aes(x = len, y = edit.distance)) + 
  geom_violin() + 
  geom_boxplot(width=0.1) +
  geom_hline(yintercept = mean(repeats$edit.distance))
ggsave(filename = "chr1_edit violin sizes more.pdf", width = 30, height = 6)
repeats$edit.distance = repeats$edit.distance - 1






