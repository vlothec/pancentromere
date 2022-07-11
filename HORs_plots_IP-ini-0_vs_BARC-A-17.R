library(stringr)
library(base)
library(msa)#
library(Biostrings)
library(seqinr)
library(ggplot2)

setwd("C:/Users/wlodz/Desktop/panC_analysis_athome/HORs_plots_IP-ini-0_vs_BARC-A-17")

do.only.similar.hors = FALSE

hors = read.csv(file = "horEdited_9852.scaffolds.bionano.final_BARC-A-17.ragtag_scaffolds.csv")

if(do.only.similar.hors)
{
  hors = hors[hors$total_variant <=0,]
}


chromosome = 1

repeatsA = read.csv(file = "repeats_summary_9852.scaffolds.bionano.final.csv", stringsAsFactors = FALSE)
repeatsA$class[is.na(repeatsA$class)] = ""
repeatsA = repeatsA[repeatsA$class == "aTha178",]
repeatsA$strand[repeatsA$strand == "+"] = "1"
repeatsA$strand[repeatsA$strand == "-"] = "2"

repeatsB = read.csv(file = "repeats_summary_BARC-A-17.ragtag_scaffolds.csv", stringsAsFactors = FALSE)
repeatsB$class[is.na(repeatsB$class)] = ""
repeatsB = repeatsB[repeatsB$class == "aTha178",]
repeatsB$strand[repeatsB$strand == "+"] = "1"
repeatsB$strand[repeatsB$strand == "-"] = "2"

repeatsTA = repeatsA[repeatsA$chromosome == unique(repeatsA$chromosome)[chromosome],]
repeatsTB = repeatsB[repeatsB$chromosome == unique(repeatsB$chromosome)[chromosome],]
repeatsT = rbind(repeatsTA, repeatsTB)

genomeA.name = "ini"
genomeB.name = "barc"
genomeA.start.coord = 14522669
genomeA.end.coord = 17705761
genomeB.start.coord = 14916339
genomeB.end.coord = 18968209

{
  print("calculate relative")
  repeatsT$HORcountA = 0
  repeatsT$HORcountB = 0
  repeatsT$HORcountShared = 0   #this is actually shared, not total
  
  for(i in 1 : nrow(hors))
  {
    print(i)
    if(hors$genomeA[i] != hors$genomeB[i])
    {
      for(j in 1 : hors$block.size.in.units[i])
      {
        repeatsT$HORcountShared[hors$start_A[i] + j - 1] = repeatsT$HORcountShared[hors$start_A[i] + j - 1] + 1
        repeatsT$HORcountShared[hors$start_B[i] + j - 1] = repeatsT$HORcountShared[hors$start_B[i] + j - 1] + 1
      }
    } else if(hors$genomeA[i] == repeatsA$assembly[1])
    {
      for(j in 1 : hors$block.size.in.units[i])
      {
        repeatsT$HORcountA[hors$start_A[i] + j - 1] = repeatsT$HORcountA[hors$start_A[i] + j - 1] + 1
        repeatsT$HORcountA[hors$start_B[i] + j - 1] = repeatsT$HORcountA[hors$start_B[i] + j - 1] + 1
      }
    } else
    {
      for(j in 1 : hors$block.size.in.units[i])
      {
        repeatsT$HORcountB[hors$start_A[i] + j - 1] = repeatsT$HORcountB[hors$start_A[i] + j - 1] + 1
        repeatsT$HORcountB[hors$start_B[i] + j - 1] = repeatsT$HORcountB[hors$start_B[i] + j - 1] + 1
      }
    }
  }
  
  repeatsT$HORrelativeSharedToTotalA = repeatsT$HORcountShared / (repeatsT$HORcountA)
  repeatsT$HORrelativeSharedToTotalB = repeatsT$HORcountShared / (repeatsT$HORcountB)
  
  repeatsT$HORrelativeSharedToTotalA[is.na(repeatsT$HORrelativeSharedToTotalA)] = 0
  repeatsT$HORrelativeSharedToTotalB[is.na(repeatsT$HORrelativeSharedToTotalB)] = 0
  repeatsT$HORrelativeSharedToTotalA[(repeatsT$HORrelativeSharedToTotalA) > 1] = 1
  repeatsT$HORrelativeSharedToTotalB[(repeatsT$HORrelativeSharedToTotalB) > 1] = 1
  
  
  png(filename = paste("Repeat_HOR_score_chr", "1", ".png", sep = ""), width = 5000, height = 3000, pointsize = 65)
  par(mfrow=c(2,2))
  smoothingSpline = smooth.spline(repeatsT$start[repeatsT$assembly == repeatsA$assembly[1]], repeatsT$HORrelativeSharedToTotalA[repeatsT$assembly == repeatsA$assembly[1]], spar = 0.3)
  plot(x = repeatsT$start[repeatsT$assembly == repeatsA$assembly[1]], 
       y = repeatsT$HORrelativeSharedToTotalA[repeatsT$assembly == repeatsA$assembly[1]],
       xlim = c(genomeA.start.coord, genomeA.end.coord),
       ylim = c(0,1),
       xlab = "Coordinates, bp",
       ylab = "repeat HOR shared to own",
       main = genomeA.name, 
       pch = ".")
  lines(smoothingSpline, col = "red", lwd = 7, lty = 3)
  
  hist(repeatsT$HORrelativeSharedToTotalA[repeatsT$assembly == repeatsA$assembly[1]],
       breaks = "Scott",
       main = "repeats HOR share histogram",
       xlab = paste(genomeA.name, " repeats HOR share", sep = ""),
       ylab = "Frequency")
  
  smoothingSpline = smooth.spline(repeatsT$start[repeatsT$assembly == repeatsB$assembly[1]], repeatsT$HORrelativeSharedToTotalB[repeatsT$assembly == repeatsB$assembly[1]], spar = 0.3)
  plot(x = repeatsT$start[repeatsT$assembly == repeatsB$assembly[1]], 
       y = repeatsT$HORrelativeSharedToTotalB[repeatsT$assembly == repeatsB$assembly[1]],
       xlim = c(genomeB.start.coord, genomeB.end.coord),
       ylim = c(0,1),
       xlab = "Coordinates, bp",
       ylab = "repeat HOR shared to own",
       main = genomeB.name, 
       pch = ".")
  lines(smoothingSpline, col = "red", lwd = 7, lty = 3)
  hist(repeatsT$HORrelativeSharedToTotalB[repeatsT$assembly == repeatsB$assembly[1]],
       breaks = "Scott",
       main = "repeats HOR share histogram",
       xlab = paste(genomeB.name, " repeats HOR share", sep = ""),
       ylab = "Frequency")
  dev.off()
}


{
  print("plot HORS")
  png(filename = paste("HOR_plot_chr", chromosome, ".png", sep = ""), width = 5000, height = 5500, pointsize = 65)
  par(mfrow=c(2,2))
  #plot1: both
  horsPlot = hors[hors$genomeA != hors$genomeB,]
  plot(xlab = genomeA.name, 
       ylab = genomeB.name, 
       x = horsPlot$start.A.bp, pch = 19, cex = 0.1, 
       y = horsPlot$start.B.bp, 
       xlim = c(genomeA.start.coord, genomeA.end.coord), 
       ylim = c(genomeB.start.coord, genomeB.end.coord),
       main = paste("HOR plot ", genomeA.name, " vs ", genomeB.name, " chr", chromosome, sep = ""))
  #plot2: genomeB
  horsPlot = hors[hors$genomeA == hors$genomeB,]
  horsPlot = horsPlot[horsPlot$genomeA == repeatsB$assembly[1],]
  plot(xlab = genomeB.name, 
       ylab = genomeB.name, 
       x = horsPlot$start.A.bp, pch = 19, cex = 0.1, 
       y = horsPlot$start.B.bp, 
       xlim = c(genomeB.start.coord, genomeB.end.coord), 
       ylim = c(genomeB.start.coord, genomeB.end.coord),
       main = paste("HOR plot ", genomeB.name, " vs ", genomeB.name, " chr", chromosome, sep = ""))
  
  #plot3: genomeA
  horsPlot = hors[hors$genomeA == hors$genomeB,]
  horsPlot = horsPlot[horsPlot$genomeA == repeatsA$assembly[1],]
  plot(xlab = genomeA.name, 
       ylab = genomeA.name, 
       x = horsPlot$start.A.bp, pch = 19, cex = 0.1, 
       y = horsPlot$start.B.bp, 
       xlim = c(genomeA.start.coord, genomeA.end.coord), 
       ylim = c(genomeA.start.coord, genomeA.end.coord),
       main = paste("HOR plot ", genomeA.name, " vs ", genomeA.name, " chr", chromosome, sep = ""))
  
  hist.hor = hist(hors$block.size.in.units, breaks = max(hors$block.size.in.units), plot = FALSE)
  try(plot(y = hist.hor$count, x = hist.hor$mids, log = "y", 
           xlab = "HOR monomer size", 
           ylab = "frequency, log scale",
           main = "HOR monomer size histogram", 
           type='h', lwd = 25, col = "pink", 
           xlim = c(0, max(hors$block.size.in.units))), silent = TRUE)
  dev.off()
}



print("plot HORS2")
png(filename = paste("BigHOR_plot_chr", chromosome, ".png", sep = ""), width = 5000, height = 5500, pointsize = 65)
par(mfrow=c(2,2))
#plot1: both
horsPlot = hors[hors$genomeA != hors$genomeB,]
plot(xlab = genomeA.name, 
     ylab = genomeB.name, 
     x = horsPlot$start.A.bp, pch = 19, cex = 1, col = "blue", 
     y = horsPlot$start.B.bp, 
     xlim = c(genomeA.start.coord, genomeA.end.coord), 
     ylim = c(genomeB.start.coord, genomeB.end.coord),
     main = paste("HOR plot ", genomeA.name, " vs ", genomeB.name, " chr", chromosome, sep = ""))
#plot2: genomeB
horsPlot = hors[hors$genomeA == hors$genomeB,]
horsPlot = horsPlot[horsPlot$genomeA == repeatsB$assembly[1],]
plot(xlab = genomeB.name, 
     ylab = genomeB.name, 
     x = horsPlot$start.A.bp, pch = 19, cex = 1, col = "blue", 
     y = horsPlot$start.B.bp, 
     xlim = c(genomeB.start.coord, genomeB.end.coord), 
     ylim = c(genomeB.start.coord, genomeB.end.coord),
     main = paste("HOR plot ", genomeB.name, " vs ", genomeB.name, " chr", chromosome, sep = ""))

#plot3: genomeA
horsPlot = hors[hors$genomeA == hors$genomeB,]
horsPlot = horsPlot[horsPlot$genomeA == repeatsA$assembly[1],]
plot(xlab = genomeA.name, 
     ylab = genomeA.name, 
     x = horsPlot$start.A.bp, pch = 19, cex = 1, col = "blue", 
     y = horsPlot$start.B.bp, 
     xlim = c(genomeA.start.coord, genomeA.end.coord), 
     ylim = c(genomeA.start.coord, genomeA.end.coord),
     main = paste("HOR plot ", genomeA.name, " vs ", genomeA.name, " chr", chromosome, sep = ""))

hist.hor = hist(hors$block.size.in.units, breaks = max(hors$block.size.in.units), plot = FALSE)
try(plot(y = hist.hor$count, x = hist.hor$mids, log = "y", 
         xlab = "HOR monomer size", 
         ylab = "frequency, log scale",
         main = "HOR monomer size histogram", 
         type='h', lwd = 25, col = "pink", 
         xlim = c(0, max(hors$block.size.in.units))), silent = TRUE)
dev.off()












pdf(file = "ini_vs_barc.pdf", width = 7.5, height = 8)
#(filename = paste("BigHOR_plot_chr", chromosome, ".png", sep = ""), width = 5000, height = 5500, pointsize = 65)
#plot1: both
horsPlot = hors[hors$genomeA != hors$genomeB,]
plot(xlab = genomeA.name, 
     ylab = genomeB.name, 
     x = horsPlot$start.A.bp, pch = 19, cex = 0.8, col = "red", 
     y = horsPlot$start.B.bp, 
     xlim = c(genomeA.start.coord, genomeA.end.coord), 
     ylim = c(genomeB.start.coord, genomeB.end.coord),
     main = paste("HOR plot ", genomeA.name, " vs ", genomeB.name, " chr", chromosome, sep = ""))
dev.off()























pdf(file = paste("Repeat_HOR_score_chr", "1", ".pdf", sep = ""), width = 30, height = 15, pointsize = 20)
par(mfrow=c(2,2))
smoothingSpline = smooth.spline(repeatsT$start[repeatsT$assembly == repeatsA$assembly[1]], repeatsT$HORrelativeSharedToTotalA[repeatsT$assembly == repeatsA$assembly[1]], spar = 0.3)
plot(x = repeatsT$start[repeatsT$assembly == repeatsA$assembly[1]], 
     y = repeatsT$HORrelativeSharedToTotalA[repeatsT$assembly == repeatsA$assembly[1]],
     xlim = c(genomeA.start.coord, genomeA.end.coord),
     ylim = c(0,1),
     xlab = "Coordinates, bp",
     ylab = "repeat HOR shared to own",
     main = genomeA.name, 
     pch = ".")
lines(smoothingSpline, col = "red", lwd = 5, lty = 3)

hist(repeatsT$HORrelativeSharedToTotalA[repeatsT$assembly == repeatsA$assembly[1]],
     breaks = "Scott",
     main = "repeats HOR share histogram",
     xlab = paste(genomeA.name, " repeats HOR share", sep = ""),
     ylab = "Frequency")

smoothingSpline = smooth.spline(repeatsT$start[repeatsT$assembly == repeatsB$assembly[1]], repeatsT$HORrelativeSharedToTotalB[repeatsT$assembly == repeatsB$assembly[1]], spar = 0.3)
plot(x = repeatsT$start[repeatsT$assembly == repeatsB$assembly[1]], 
     y = repeatsT$HORrelativeSharedToTotalB[repeatsT$assembly == repeatsB$assembly[1]],
     xlim = c(genomeB.start.coord, genomeB.end.coord),
     ylim = c(0,1),
     xlab = "Coordinates, bp",
     ylab = "repeat HOR shared to own",
     main = genomeB.name, 
     pch = ".")
lines(smoothingSpline, col = "red", lwd = 5, lty = 3)
hist(repeatsT$HORrelativeSharedToTotalB[repeatsT$assembly == repeatsB$assembly[1]],
     breaks = "Scott",
     main = "repeats HOR share histogram",
     xlab = paste(genomeB.name, " repeats HOR share", sep = ""),
     ylab = "Frequency")
dev.off()


