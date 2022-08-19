
#!/usr/bin/env Rscript
thisFile = function() 
{
  cmd.Args = commandArgs(trailingOnly = FALSE)
  find.file = "--file="
  match = grep(find.file, cmd.Args)
  if (length(match) > 0) {
    return(normalizePath(sub(find.file, "", cmd.Args[match])))
  } else {
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

installation.path = thisFile()
installation.path = strsplit(installation.path, split = "/")[[1]]
installation.path = paste(installation.path[1:(length(installation.path) - 2)], collapse = "/")
execution.path = getwd()

lib.path = paste(installation.path, "/libs", sep = "")
#installation.path = "/rds/user/pw457/hpc-work/TRASH"
#execution.path = getwd()
#lib.path = "/rds/user/pw457/hpc-work/TRASH/libs"
#fasta.list = c("51.atTest1.fa", "50.atTest1.fasta")
.libPaths(lib.path)


library(stringr, quietly = TRUE)
library(base, quietly = TRUE)
library(Biostrings, quietly = TRUE)
library(seqinr, quietly = TRUE)
library(doParallel, quietly = TRUE)

src.files = list.files(path = paste(installation.path, "/src", sep = ""), pattern = "fn_", full.names = TRUE)

for(i in 1 : length(src.files))
{
  print(paste("attaching function from: ", src.files[i], sep = ""))
  source(src.files[i])
}
arguments = commandArgs(trailingOnly = TRUE)

#run settings
{
  sequence.templates = NA
  skip.repetitive.regions = FALSE
  change.lib.paths = TRUE
  cores.no = NA
  set.kmer = 12 # kmer size used for initial identification of repetitive regions
  set.threshold = 10 # window repetitiveness score (0-100) threshold
  set.max.repeat.size = 1000 # max size of repeats to be identified
  filter.small.regions = 1500 # repetitive windows smaller than this size will be removed (helps getting rid of regions with short duplications)
  filter.small.repeats = 4 # repetitive windows where dominant kmer distance is lower than this value will be removed (for example AT dinucleotide repeats)
  window.size = 1500 # how far apart kmers can be in the initial search for exact matches. No repeats larger than this will be identified
}
fasta.list = NULL
{
  if(length(arguments) == 0)
  {
    stop("At least one argument must be supplied", call. = FALSE)
  } else
  {
    for(i in 1 : length(arguments))
    {
      if(arguments[i] == "--def")
      {
        change.lib.paths = FALSE
      } else if(arguments[i] == "--skipr")
      {
        skip.repetitive.regions = TRUE
      } else if(arguments[i] == "-w")
      {
        set.kmer = as.numeric(arguments[i + 1])
      } else if(arguments[i] == "-t")
      {
        set.threshold = as.numeric(arguments[i + 1])
      } else if(arguments[i] == "-m")
      {
        set.max.repeat.size = as.numeric(arguments[i + 1])
      } else if(arguments[i] == "-freg")
      {
        filter.small.regions = as.numeric(arguments[i + 1])
      } else if(arguments[i] == "-frep")
      {
        filter.small.repeats = as.numeric(arguments[i + 1])
      } else if(arguments[i] == "-o")
      {
        execution.path = arguments[i + 1]
      } else if(arguments[i] == "-win")
      {
        window.size = as.numeric(arguments[i + 1])
      }
      else if(arguments[i] == "-seqt")
      {
        sequence.templates = arguments[i + 1]
      }
      else if(arguments[i] == "-par")
      {
        cores.no = as.numeric(arguments[i + 1])
      }
      if(grepl(".fa", arguments[i]))
      {
        if(file.exists(paste(execution.path, arguments[i], sep = "/")))
        {
          fasta.list = c(paste(execution.path, arguments[i], sep = "/"))
        } else if(file.exists(arguments[i]))
        {
          fasta.list = c(fasta.list, arguments[i])
        } else
        {
          stop("Cannot find the", arguments[i], "fasta file specified", call. = FALSE)
        }
      }
    }
  }
  if(Sys.info()['sysname'] != "Linux")
  {
    print("Only Linux OS tested")
  } 
}

# read in fasta and setup parallel to the length of fasta file
print(fasta.list)
{
  fasta.sequence = NULL
  sequences = data.frame(file.name = "", fasta.name = "")
  for(i in 1 : length(fasta.list))
  {
    fasta = read.fasta(file = fasta.list[i], seqtype = "DNA", as.string = TRUE)
    for(j in 1 : length(fasta))
    {
      fasta.sequence = c(fasta.sequence,toupper(fasta[j]))
      sequences = rbind(sequences, data.frame(file.name = strsplit(fasta.list[i], split = "/")[[1]][length(strsplit(fasta.list[i], split = "/")[[1]])], fasta.name = names(fasta[j])))
      strsplit(fasta.list[i], split = "/")[[1]]
    }
  }
  sequences = sequences[-1,]
  if(!is.na(cores.no))
  {
    registerDoParallel(cores = cores.no)
  } else
  {
    registerDoParallel(cores = nrow(sequences))
  }
  print("Currently registered parallel backend name, version and cores")
  print(getDoParName())
  print(getDoParVersion())
  print(getDoParWorkers())
}

if(!is.na(sequence.templates)) 
{
  sequence.templates = read.csv(file = sequence.templates)
}

#identify repeats per chromosome
foreach(i = 1 : length(fasta.sequence)) %dopar% {
  #for(i in 1 : length(fasta.sequence)) {
  print(paste("Working on sequence ", i , sep = ""))
  if(!skip.repetitive.regions)
  {
    repetitive.regions = Repeat.Identifier(DNA.sequence = fasta.sequence[i], assemblyName = sequences$file.name[i], fasta.name = sequences$fasta.name[i], 
                                           kmer = set.kmer, window = window.size, threshold = set.threshold, mask.small.regions = filter.small.regions, mask.small.repeats = filter.small.repeats,
                                           max.repeat.size = set.max.repeat.size,
                                           tests = 4, temp.folder = execution.path, sequence.template = sequence.templates, mafft.bat.file = paste(installation.path, "/src/mafft-linux64/mafft.bat", sep = ""))
    write.csv(x = repetitive.regions, file = paste(execution.path, "/", sequences$file.name[i], "_out/Regions_", sequences$file.name[i], "_", sequences$fasta.name[i], ".csv", sep = ""), row.names = FALSE)
  }
  gc()
} 

#save repeats per genome
for(i in 1 : length(fasta.list))
{
  setwd(paste(execution.path, "/", sep = ""))
  repeat.files = system(paste("find . -name \"Regions*\"", sep = ""), intern = TRUE)
  regions = NULL
  for(j in 1 : length(repeat.files))
  {
    region.read = read.csv(file = paste(repeat.files[j], sep = ""))
    if(nrow(region.read) > 1)
    {
      regions = rbind(regions, region.read)
    }
  }
  regions = regions[,-c(1)]
  write.csv(x = regions, file = paste(execution.path, "/Summary.of.repetitive.regions.", sequences$file.name[i], ".csv", sep = ""), quote = FALSE)
  extract.all.repeats(temp.folder = execution.path, 
                      assemblyName = strsplit(fasta.list[i], split = "/")[[1]][length(strsplit(fasta.list[i], split = "/")[[1]])])
  
  draw.scaffold.repeat.plots(temp.folder = execution.path, 
                             assemblyName = strsplit(fasta.list[i], split = "/")[[1]][length(strsplit(fasta.list[i], split = "/")[[1]])], 
                             fastaDirectory = fasta.list[i], 
                             only.Chr1_5 = FALSE, 
                             single.pngs = TRUE)
}

#do HORs per chromosome
foreach(i = 1 : length(fasta.sequence)) %dopar% {
  
  #extract all repeats that have a class assigned
  
  #HOR.wrapper(temp.folder = execution.path, assemblyName = sequences$file.name[i], mafft.bat.file = paste(installation.path, "/src/mafft-linux64/mafft.bat", sep = ""))
  #gc()
} 



{
  #sequence.templates = "~/sequence.template.csv" # path to a csv file used to match identified repeats to templates so they are globally in the same frame (same start position), set as NA to skip
  
  
  #if(!is.na(sequence.templates)) # modify which templates should be used
  #{
  #  sequence.templates = read.csv(file = sequence.templates)
  #  sequence.templates = sequence.templates[sequence.templates$group == "ath",]
  #}
  
}
