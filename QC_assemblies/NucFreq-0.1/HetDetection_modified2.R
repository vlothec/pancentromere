#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(directlabels)
library(plyr)
require(gridExtra)
require(scales)
library(RColorBrewer)
library(reshape2)
library(data.table)
library(dplyr)
library(tidyr)
library(zoo)

###########################################################################
###########################################################################
## The original script 'HetDetection.R' has been modified in this study  ##
## to make it executable from command line, and to fix the incorrect     ##
## merging of blocks across different chromosomes                        ##
###########################################################################
###########################################################################

#clear variables
#rm(list=ls(all=TRUE))

#set working directory
setwd(args[1])

#load NucFreq bed file
df = fread(args[2], stringsAsFactors = FALSE, fill=TRUE, quote="", header=FALSE, skip=2)
cols =  c("chr", "start", "end", "first", "second")
colnames(df) <- cols

#determine the ratio of the first and second most common bases
df$het_ratio = round(df$second/(df$first+df$second)*100, 1)

#filter if the het ratio is >= 10%
df1 = df %>% 
  group_by(chr) %>% 
  filter(het_ratio >= 10)

#calculate the distance (in bp) between consecutive positions
df2 = df1 %>%
  group_by(chr) %>%
  mutate(distance = start - lag(start, default = start[1]))

#shift the distance column up one row
shift <- function(x, n){
  c(x[-(seq(n))], rep(NA, n))
}

df2$distance <- shift(df2$distance, 1)

#filter rows with a distance <=500 bp between positions (i.e. the het must have another base change within 500 bp)
df3 = df2 %>% 
  group_by(chr) %>% 
  filter(distance <= 500)
df3 = df3 %>%
  group_by(chr) %>%
  mutate(distance2 = start - lag(start, default = start[1]))
df3$distance2 <- shift(df3$distance2, 1)

#duplicate top row and change its distance to 501 bp (to get rows in register)
df4 = df3 %>% 
  group_by(chr) %>% 
  filter(row_number() <= 1) %>% 
  bind_rows(df3)
df5 = df4 %>%
  arrange(start, .by_group = TRUE) %>%
  mutate(distance2 = replace(distance2, row_number() == 1, 501))

#shift up the end column to get the range of the hets on one row
df5$end <- shift(df5$end, 1)

#filter only if there are 5 consecutive rows of distance <=500 bp (i.e. the het must have 5 base changes within 500 bp)
r <- with(with(df5, rle(distance2<=500)),rep(lengths,lengths))
df5$het <- with(df5,distance2<=500) & (r>=4)

###### Modifications implemented by Fernando Rabanal start here ########
###### Added 12.07.2022 ########

# Define blocks of heterozygosity
tmp_block <- 1
blocks <- 1
for(i in 1:(length(df5$het)-1)){
  if(is.na(df5$het[i]==df5$het[i+1])){
    break
  }else if(df5$het[i]==df5$het[i+1]){
    tmp_block <- tmp_block
  }else{
    tmp_block <- i+1
  }
  blocks <- c(blocks, tmp_block)
}
df5$blocks <- blocks

# Filter rows defined as het
df6 <- filter(df5, het == "TRUE")

# Filter out distance2 == 0 because their coordinates come from different chromosomes 
df6 = df6 %>% 
  group_by(chr) %>% 
  filter(distance2 > 0)

# 
block_ids <- unique(df6$blocks)

HET_df <- data.frame(chr=NA, start=NA, end=NA, het_ratio=NA, het_length=NA)

if(dim(df6)[1] > 0){
  for(i in 1:length(block_ids)){
    df_tmp <- subset(df6, blocks == block_ids[i])
  
    HET_df[i,1] <- df_tmp[1,1]
    HET_df[i,2] <- df_tmp[1,2]
    HET_df[i,3] <- df_tmp[dim(df_tmp)[1],3]
    HET_df[i,4] <- round(mean(df_tmp$het_ratio), digits=1)
    HET_df[i,5] <- HET_df[i,3] - HET_df[i,2]
  }  
}else{
  HET_df <- HET_df
}

#print the table
write.table(HET_df, paste(args[2], ".het_regions.tbl", sep=""), row.names = F, quote = F, sep="\t")


