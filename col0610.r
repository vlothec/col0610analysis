#
# R 4.0.3 
# updated 22 Feb 2021
# 

#####
# 0. Functions and libraries
revCompString = function(DNAstr) {  return(toupper(toString(reverseComplement(DNAString(DNAstr))))) }
library(stringr)
library(beepr)
library(svMisc)
library(base)
library(magrittr)
library(msa)
library(Biostrings)
library(seqinr)

remove(frame1, frame2, frame3, frame4, frame5, frame6, genome)

genome = "C:/Users/wlodz/Desktop/t2t-col.20210610.fasta"

time1 = Sys.time()
frame1 = identify.repetitive.regions(fastaDirectory = genome, kmer = 10, window = 1000, threshold = 10, filter.small = 3500)
a = Sys.time() - time1

time2 = Sys.time()
frame2 = closest.identical.kmer.distances(regions.data.frame = frame1, kmer = 6, reach = 300, filter.smaller.N = 4, plot = TRUE)
b = Sys.time() - time2


time3 = Sys.time()
frame3 = test.random.Nmers(regions.data.frame = frame2, tests = 3)
c = Sys.time() - time3



time4 = Sys.time()
revCompString = function(DNAstr) {  return(toupper(toString(reverseComplement(DNAString(DNAstr))))) }
frame4 = generate.secondary.consensus(regions.data.frame = frame3)
d = Sys.time() - time4

frame7 = frame4[,-c(6)]

write.csv(x = frame7, file = "frame7col0610.csv", quote = FALSE)

a
b
c
d

time5 = Sys.time()
frame5 = merge.and.generate.shifted.consensus(regions.data.frame = frame4, merge.threshold = 0.65)
e = Sys.time() - time5


e

frame5dv = frame5

time6 = Sys.time()
frame6 = match.repeats(regions.data.frame = frame5, fastaDirectory = genome)
f = Sys.time() - time6


f



#####
# 1. Identify repetitive regions

identify.repetitive.regions = function(fastaDirectory = "C:/Users/wlodz/Desktop/testmanuscript/test.fasta", kmer = 12, window = 1000, threshold = 10, filter.small = 0)
{
  time = Sys.time()
  
  fasta = read.fasta(fastaDirectory, as.string = TRUE)
  
  repetitive.windows = data.frame(index = integer(), fasta.name = character(), start = integer(), end = integer(), average.score = double(), window.sequence = character(), most.freq.value.N = integer())
  
  index = 0
  
  for(i in 1 : length(fasta))   # for every fasta sequence in the file
  {
    
    seqA = fasta[[i]]
    seqA = toupper(seqA)
    seq.length = nchar(seqA)
    scores = vector(mode = "numeric", length = seq.length %/% window + 1)
    
    for(ii in 1 : (length(scores) - 1))  # for every window of size window but not the last one
    {
      print(paste("sequence ", i, "/", length(fasta), ", window ", ii, "/", length(scores), sep = ""))
      window.start = (ii - 1) * window
      check.window = str_sub(seqA, window.start, (window.start + window - 1))
      for(iii in 1 : (window - kmer))
      {
        check.string = str_sub(string = check.window, start = iii, end = (iii + kmer - 1))
        count = 1 * ((str_count(string = check.string, pattern = "N") == 0) & str_count(string = check.window, pattern = check.string) > 1)
        scores[ii] = scores[ii] + count
      }
    }
    ii = length(scores) # for the last window
    print(paste("sequence ", i, "/", length(fasta), ", window ", ii, "/", length(scores), sep = ""))
    window.start = (ii - 1) * window
    check.window = str_sub(seqA, window.start)
    for(iii in 1 : (nchar(check.window) - kmer))
    {
      check.string = str_sub(string = check.window, start = iii, end = (iii + kmer - 1))
      count = 1 * ((str_count(string = check.string, pattern = "N") == 0) & str_count(string = check.window, pattern = check.string) > 1)
      scores[ii] = scores[ii] + count
    }
    
    scores[1 : (length(scores) - 1)] = scores[1 : (length(scores) - 1)] / window * 100 #change score into a percentage score
    scores[length(scores)] = scores[length(scores)] / nchar(check.window)
    
    ii = 1
    while(ii <= length(scores)) #apply threshold and save repetitive windows from this fasta
    {
      if(scores[ii] >= threshold)
      {
        index = index + 1
        name = names(fasta)[i]
        start = (ii - 1) * window
        end = start - 1
        ave.score = 0
        s = 0
        while(scores[ii] >= threshold)
        {
          s = s + 1
          end = end + 1000
          ave.score = scores[ii] + ave.score 
          ii = ii + 1
        }
        ave.score = ave.score / s
        window.sequence = str_sub(seqA, start, end)
        most.freq.value.N = 0
        temp = data.frame(index, name, start, end, ave.score, window.sequence, most.freq.value.N)
        repetitive.windows = rbind(repetitive.windows, temp)
      }
      ii = ii + 1
    }
    
    repetitive.windows = repetitive.windows[(repetitive.windows$end - repetitive.windows$start) > filter.small,]
    
    print(paste("all done in ", (Sys.time() - time), sep = ""))
  }
  gc()
  return(repetitive.windows)
}





#####
# 2. Calculate closest identical k-mer distances

closest.identical.kmer.distances = function(regions.data.frame, kmer = 12, reach = 5000, filter.smaller.N = 0, plot = TRUE)
{
  for(i in 1 : nrow(regions.data.frame))
  {
    seqB = regions.data.frame$window.sequence[i]
    print(paste("calculating distances on region ", i, "/", nrow(regions.data.frame), " region size: ", nchar(seqB), sep = ""))
    distance = vector(mode = "numeric", length =  nchar(seqB))
    for(ii in 1 : (nchar(seqB) - kmer + 1))
    {
      kmer.pattern = str_sub(seqB, ii, ii + kmer - 1)
      window.string = str_sub(seqB, ii + 1, ii + reach)
      distance[ii] = str_locate(string = window.string, pattern = kmer.pattern)[[1]]
    }
    distance = distance[!is.na(distance)]
    distance = distance[distance > 0]
    if(length(distance) > 0)
    {
      a = hist(distance, breaks = max(distance), plot = FALSE)
      regions.data.frame$most.freq.value.N[i] = a$breaks[which.max(a$counts) + 1]
      if(plot == TRUE)
      {
        png(filename = paste("C:/Users/wlodz/Desktop/histograms/histogram.of.", regions.data.frame$name[i], ".", regions.data.frame$index[i], ".N=", regions.data.frame$most.freq.value.N[i], ".png", sep = ""), width = 750, height = 500, pointsize = 25)
        hist(distance, breaks = max(distance), main = paste("Histogram of", regions.data.frame$name[i], "index", regions.data.frame$index[i], "N =", regions.data.frame$most.freq.value.N[i], sep = " "))
        dev.off()
      }
    }
  }
  regions.data.frame = regions.data.frame[regions.data.frame$most.freq.value.N >= filter.smaller.N,]
  gc()
  return(regions.data.frame)
}




#####
# 3. Test random N-size mers, handle overlaps, choose the best one, extract repeats for mafft alignment

test.random.Nmers = function(regions.data.frame, tests = 5)
{
  consensus.primary = vector(mode = "character", length = nrow(regions.data.frame))
  consensus.count = vector(mode = "numeric", length = nrow(regions.data.frame))
  for(i in 1 : nrow(regions.data.frame))
  {
    print(paste("testing region ", i, "/", nrow(regions.data.frame), " window size: ", regions.data.frame$end[i] - regions.data.frame$start[i], sep = ""))
    seqC = regions.data.frame$window.sequence[i]
    N = regions.data.frame$most.freq.value.N[i]
    if(N == 0)
    {
      N = 1
    }
    random.sequences.start = sample(1 : (nchar(seqC) - N), tests)
    random.sequences = str_sub(seqC, random.sequences.start, (random.sequences.start + N - 1))
    random.sequence.scores = vector(mode = "numeric", length = tests)
    random.sequence.counts = vector(mode = "numeric", length = tests)
    
    maxMis = N %/% 3
    if(maxMis > 100)
    {
      maxMis = 100
    }
    
    #match = countPDict(pdict = random.sequences, subject = seqC, max.mismatch = maxMis, with.indels = TRUE)#check
    
    for(ii in 1 : tests)
    {
      match = matchPattern(pattern = random.sequences[ii], subject = seqC, max.mismatch = maxMis, with.indels = TRUE)
      temp = as.data.frame(match)
      temp$start  = start(match)
      temp$end = end(match)
      temp$strand = "+"
      iii = nrow(temp)
      random.sequence.counts[ii] = 0
      if(iii > 0)
      {
        random.sequence.counts[ii] = nrow(temp) 
      }
      while(iii > 1)
      {
        if(temp$start[iii] - temp$end[iii - 1] < 0)
        {
          temp = temp[-iii,]
        }
        iii = iii - 1
      }
      random.sequence.scores[ii] = sum(temp$end - temp$start)
      
    }
    consensus.primary[i] = random.sequences[which.max(random.sequence.scores)]
    consensus.count[i] = max(random.sequence.counts)
  }
  gc()
  return(cbind(regions.data.frame, consensus.primary, consensus.count))
}


#####
# 4. Refine each consensus by mapping again, aligning with mafft and extracting consensus 

generate.secondary.consensus = function(regions.data.frame, assemblyName = "")
{
  if(assemblyName == "")
  {
    assemblyName = date()
  }
  consensus.secondary = vector(mode = "character", length = nrow(regions.data.frame))
  repeats.identified = vector(mode = "numeric", length = nrow(regions.data.frame))
  
  for(i in 1 : nrow(regions.data.frame))
  {
    print(paste("generating consensus for region ", i, "/", nrow(regions.data.frame), " window size: ", regions.data.frame$end[i] - regions.data.frame$start[i], sep = ""))
    seqC = regions.data.frame$window.sequence[i]
    N = regions.data.frame$most.freq.value.N[i]
    maxMis = N %/% 3
    if(maxMis > 100)
    {
      maxMis = 100
    }
    matchPlus = matchPattern(pattern = regions.data.frame$consensus.primary[i], subject = seqC, max.mismatch = maxMis, with.indels = TRUE)
    tempP = as.data.frame(matchPlus)
    tempP$start  = start(matchPlus)
    tempP$end = end(matchPlus)
    if(nrow(tempP) > 0)
    {
      tempP$strand = "+"
    }
    matchMinus = matchPattern(pattern = toupper(revCompString(regions.data.frame$consensus.primary[i])), subject = seqC, max.mismatch = maxMis, with.indels = TRUE)
    tempM = as.data.frame(matchMinus)
    tempM$start  = start(matchMinus)
    tempM$end = end(matchMinus)
    if(nrow(tempM) > 0)
    {
      tempM$strand = "-"
    }
    match = rbind(tempP, tempM)
    match = match[order(match$start),]
    
    ##########
    iii = nrow(match)
    
    while(iii > 1)
    {
      overlap = match$start[iii] - match$end[iii - 1] - 1
      if(overlap < 0 & overlap >= -10)
      {
        overlap = -overlap
        b = overlap %/% 2
        c = overlap - b
        match$start[iii] = match$start[iii] + b
        match$end[iii - 1] = match$end[iii - 1] - c
        
      } else if(overlap < -10 & overlap > -500)
      {
        match = match[-c(iii),] 
      }
      iii = iii - 1
    }
    
    setwd("C:/Users/wlodz/Desktop/mafft-7.475-win64-signed/mafft-win")
    dir.create(paste("C:/Users/wlodz/Desktop/mafft-7.475-win64-signed/mafft-win/", assemblyName, sep = ""))
    dir.create(paste("C:/Users/wlodz/Desktop/mafft-7.475-win64-signed/mafft-win/", assemblyName, "/inputs", sep = ""))
    dir.create(paste("C:/Users/wlodz/Desktop/mafft-7.475-win64-signed/mafft-win/", assemblyName, "/outputs", sep = ""))
    
    if(nrow(match) > 0)
    {
      write.csv(x = match, file = paste(assemblyName, "/Repeats_", regions.data.frame$name[i], "_", regions.data.frame$start[i], "_", regions.data.frame$end[i], ".csv", sep = ""))
      #export sequences in fasta
      for(ii in 1 : nrow(match))
      {
        if(match$strand[ii] == "+")
        {
          write.fasta(sequences = match$x[ii], names = paste("primary.extract", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".", ii, sep = ""), file.out = paste("C:/Users/wlodz/Desktop/mafft-7.475-win64-signed/mafft-win/", assemblyName, "/inputs/primary.extract", ".", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".", "fasta", sep = ""), open = "a", as.string = TRUE)
        } else if(match$strand[ii] == "-")
        {
          write.fasta(sequences = revCompString(match$x[ii]), names = paste("primary.extract", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".", ii, sep = ""), file.out = paste("C:/Users/wlodz/Desktop/mafft-7.475-win64-signed/mafft-win/", assemblyName, "/inputs/primary.extract", ".", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".", "fasta", sep = ""), open = "a", as.string = TRUE)
        }
      }
      
      #mafft
      
      input = paste("C:/Users/wlodz/Desktop/mafft-7.475-win64-signed/mafft-win/", assemblyName, "/inputs/primary.extract", ".", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".", "fasta", sep = "")
      output = paste("C:/Users/wlodz/Desktop/mafft-7.475-win64-signed/mafft-win/", assemblyName, "/outputs/primary.extract", ".", regions.data.frame$index[i], ".", regions.data.frame$name[i], ".", "aligned.fasta", sep = "")
      system(paste("mafft.bat --retree 2 --inputorder ", input, " > ", output, sep = ""))
      
      alignment = read.alignment(output, format = "FASTA", forceToLower = FALSE)
      consensus = consensus(alignment, threshold = 0.6)
      consensus = toupper(consensus[consensus != "-"])
      consensus = paste(consensus, collapse = "")
      repeats.identified[i] = length(alignment$seq)
      
      consensus.secondary[i] = consensus
      remove(consensus)
    }
    
  }
  gc()
  return(cbind(regions.data.frame, consensus.secondary, repeats.identified))
}


#####
# 5. Align secondary consensus sequences to classify and shift, replace with merged as tertiary consensus

merge.and.generate.shifted.consensus = function(regions.data.frame, merge.threshold = 0.65)
{
  consensus.shifted = vector(mode = "character", length = nrow(regions.data.frame))
  consensus.second.doubled = vector(mode = "character", length = nrow(regions.data.frame))
  
  regions.data.frame$class = 0
  
  class = 1
  
  for(i in 1 : nrow(regions.data.frame))
  {
    #print(i)
    consensus.shifted[i] = regions.data.frame$consensus.secondary[i]
    consensus.second.doubled[i] = paste(regions.data.frame$consensus.secondary[i], regions.data.frame$consensus.secondary[i], sep = "")
  }
  
  #matrix.scores = matrix(nrow = nrow(regions.data.frame), ncol = nrow(regions.data.frame))
  
  for(i in 1 : (nrow(regions.data.frame) - 1))
  #for(i in 58 : 60)
  {
    ref = pairwiseAlignment(pattern = consensus.second.doubled[i], subject = consensus.second.doubled[i] ,type = "local", scoreOnly = TRUE)
    #matrix.scores[i,i] = ref
    print(paste(i, i, ref, sep = " "))
    #for(ii in (i + 1) : 61)
    for(ii in (i + 1) : nrow(regions.data.frame))
    {
      #print(ii)
      check = pairwiseAlignment(pattern = consensus.second.doubled[i], subject = consensus.second.doubled[ii] ,type = "local", scoreOnly = TRUE)
      check2 = pairwiseAlignment(pattern = consensus.second.doubled[ii], subject = consensus.second.doubled[ii] ,type = "local", scoreOnly = TRUE)
      #matrix.scores[i,ii] = check
      print(paste(i, ii, check, sep = " "))
      if(check/ref > merge.threshold & check/check2 > merge.threshold)
      {
        if(regions.data.frame$class[i] == "0")
        {
          regions.data.frame$class[i] = class
          class = class + 1
        }
        if(regions.data.frame$class[i] != "0")
        {
          regions.data.frame$class[ii] = regions.data.frame$class[i]
        }
      }
    }
  }
  
  setwd("C:/Users/wlodz/Desktop/mafft-7.475-win64-signed/mafft-win")
  i = 1 #extract same class repeats
  while(i <= max(regions.data.frame$class))
  {
    which(regions.data.frame$class == i)
    for(ii in 1 : length(which(regions.data.frame$class == i)))
    {
      write.fasta(sequences = consensus.second.doubled[which(regions.data.frame$class == i)][ii], names = ii, file.out = paste("C:/Users/wlodz/Desktop/mafft-7.475-win64-signed/mafft-win/shifting/", i, "wip.fasta", sep = ""), open = "a", as.string = TRUE)
    }
    input = paste("C:/Users/wlodz/Desktop/mafft-7.475-win64-signed/mafft-win/shifting/", i, "wip.fasta", sep = "")
    output = paste("C:/Users/wlodz/Desktop/mafft-7.475-win64-signed/mafft-win/shifting/", i, "wip.aligned.fasta", sep = "")
    system(paste("mafft.bat --retree 2 --inputorder ", input, " > ", output, sep = ""))
    alignment = read.alignment(output, format = "FASTA", forceToLower = FALSE)

    shifts = vector(mode = "numeric", length = length(which(regions.data.frame$class == i)))
  
    for(ii in 1 : length(which(regions.data.frame$class == i)))
    {
      ali.temp = alignment$seq[ii]
      positions = c(str_locate(string = ali.temp,pattern = "a")[1], str_locate(string = ali.temp,pattern = "t")[1], str_locate(string = ali.temp,pattern = "g")[1], str_locate(string = ali.temp,pattern = "c")[1])
      shifts[ii] = min(positions)
    }
    
    for(ii in 1 : length(which(regions.data.frame$class == i)))
    {
      shift = max(shifts) - shifts[ii] + 1
      consensus.shifted[which(regions.data.frame$class == i)[ii]] = str_sub(string = consensus.second.doubled[which(regions.data.frame$class == i)[ii]], start = shift, end = (shift -1 + regions.data.frame$most.freq.value.N[which(regions.data.frame$class == i)[ii]]))
    }
    i = i + 1
  }
  for(i in 1 : nrow(regions.data.frame))
  {
    if(is.na(consensus.shifted[i]))
    {
      consensus.shifted[i] = regions.data.frame$consensus.secondary[i]
    }
  }
  gc()
  return(cbind(regions.data.frame, consensus.shifted))
}


#####
# 6. Generate new data frame with repeats classes, matchPattern them to the input sequences

match.repeats = function(regions.data.frame, fastaDirectory = "C:/Users/wlodz/Desktop/testmanuscript/test.fasta")
{
  fasta = read.fasta(fastaDirectory, as.string = TRUE)
  
  repeats.classes = data.frame(index = integer(), repeat.size = character(), number.of.repeats = integer(), average.pairwise.distance = double(), average.pairwise.distance.sd = double(), consensus.sequence = character())
  
  for(i in 1 : nrow(regions.data.frame))
  {
    print(paste("Analysing repeat class from region ", i, "/", nrow(regions.data.frame), sep = ""))

    maxMis = nchar(regions.data.frame$consensus.shifted[i]) %/% 3
    if(maxMis > 100)
    {
      maxMis = 100
    }
    seqF = regions.data.frame$window.sequence[i]
     
    matchPlus = matchPattern(pattern = regions.data.frame$consensus.shifted[i], subject = seqF, max.mismatch = maxMis, with.indels = TRUE)
    tempP = as.data.frame(matchPlus)
    tempP$start  = start(matchPlus)
    tempP$end = end(matchPlus)
    if(nrow(tempP) > 0)
    {
      tempP$strand = "+"
    }
    matchMinus = matchPattern(pattern = toupper(revCompString(regions.data.frame$consensus.shifted[i])), subject = seqF, max.mismatch = maxMis, with.indels = TRUE)
    tempM = as.data.frame(matchMinus)
    tempM$start  = start(matchMinus)
    tempM$end = end(matchMinus)
    if(nrow(tempM) > 0)
    {
      tempM$strand = "-"
    }
    match = rbind(tempP, tempM)
    match = match[order(match$start),]
    
    iii = nrow(match)  
    
    while(iii > 1)
    {
      overlap = match$start[iii] - match$end[iii - 1] - 1
      if(overlap < 0 & overlap >= -10)
      {
        overlap = -overlap
        b = overlap %/% 2
        c = overlap - b
        match$start[iii] = match$start[iii] + b
        match$end[iii - 1] = match$end[iii - 1] - c
        
      } else if(a < -10 & a > -500)
      {
        match = match[-c(i),] 
      }
      iii = iii - 1
    }
    
      
    if(length(match$x != 0))#######number of repeats
    {
      for(iii in 1 : length(match$x)) #write fasta to align
      {
        if(match$strand[iii] == "+")
        {
          write.fasta(sequences = match$x[iii], names = paste(regions.data.frame$name[i], regions.data.frame$class[i], match$start[iii], match$end[iii], match$strand[iii], "fasta", sep = "."), file.out = paste("C:/Users/wlodz/Desktop/mafft-7.475-win64-signed/mafft-win/inputs/final.extract", i, regions.data.frame$class[i], "fasta", sep = "."), open = "a", as.string = TRUE)
          
        } else if(match$strand[iii] == "-")
        {
          write.fasta(sequences = revCompString(match$x[iii]), names = paste(regions.data.frame$name[i], regions.data.frame$class[i], match$start[iii], match$end[iii], match$strand[iii], "fasta", sep = "."), file.out = paste("C:/Users/wlodz/Desktop/mafft-7.475-win64-signed/mafft-win/inputs/final.extract", i, regions.data.frame$class[i], "fasta", sep = "."), open = "a", as.string = TRUE)
        }
      }
    }
    
    
    setwd("C:/Users/wlodz/Desktop/mafft-7.475-win64-signed/mafft-win")
    input = paste("C:/Users/wlodz/Desktop/mafft-7.475-win64-signed/mafft-win/inputs/final.extract", i, regions.data.frame$class[i], "fasta", sep = ".")
    output =  paste("C:/Users/wlodz/Desktop/mafft-7.475-win64-signed/mafft-win/outputs/final.extract", i, regions.data.frame$class[i], ".aligned.fasta", sep = ".")
    system(paste("mafft.bat --sparsescore 1000 --inputorder ", input, " > ", output, sep = "")) 
      
    alignment = read.alignment(output, format = "FASTA", forceToLower = FALSE)
    consensus = consensus(alignment, threshold = 0.6)
    consensus = toupper(consensus[consensus != "-"])
    consensus = paste(consensus, collapse = "") #######consensus.sequence
    nchar(consensus) ########repeat.size
    ave.score = ave(as.matrix(dist.alignment(alignment, matrix = "identity")))[1,1] * 100   ########average.pairwise.distance
    ave.score.sd = sd(as.matrix(dist.alignment(alignment, matrix = "identity"))) * 100   ########average.pairwise.distance.sd
      
    temp = data.frame(index = regions.data.frame$class[i], repeat.size = nchar(consensus) , number.of.repeats = nrow(match), average.pairwise.distance = ave.score, average.pairwise.distance.sd = ave.score.sd, consensus.sequence = consensus)
    repeats.classes = rbind(repeats.classes, temp)
  
  }#end iterating over regions
  return(repeats.classes)
}#end function

#####
#7. matchpattern chosen repeat to the whole sequence and extract 

############


cen180plus = tolower("AGTATAAGAACTTAAACCGCAACCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAAGACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAGAAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTG")
cen180minus = tolower("CAATCATACACATGACATCAAGTCATATTCGACTCCAAAACACTAACCAACCTTCTTCTTGCTTCTCAAAGCTTTCATGGTGTAGCCAAAGTCCATATGAGTCTTTGGCTTTGTGTCTTCTAACAAGGAAACACTACTTAGGCTTTTAAGATCGGTTGCGGTTTAAGTTCTTATACT")

class.names = c("cen180plus", "cen180minus")

class.sequences = c(cen180plus, cen180minus)

allRepeats = data.frame()
allRepeats$contig = vector(mode = "numeric")
allRepeats$class = vector(mode = "numeric")
allRepeats$start = vector(mode = "numeric")
allRepeats$end = vector(mode = "numeric")
allRepeats$sequence = vector(mode = "character")
allRepeats$strand = vector(mode = "character")

a = read.fasta(file = "C:/Users/wlodz/Desktop/cen1-5.col0601.fasta", as.string = TRUE)


for(i in 1 : 5)
{
  time = Sys.time()
  print(paste("fasta file ", i, "/7", sep = ""))
  for(j in 1 : 2)#length(class.names))
  {
    print(paste("repeat class: ", class.names[j], sep = ""))
    maxMis = nchar(class.sequences[1]) %/% 3
    match = matchPattern(pattern = class.sequences[j], subject = a[[i]][1], max.mismatch = maxMis, with.indels = TRUE)
    temp = as.data.frame(match)
    if(length(temp$x) > 0)
    {
      print(i)
      print(j)
      temp$contig = getAnnot(a[[i]])
      temp$start  = start(match)
      temp$end = end(match)
      temp$strand = '+'
      temp$class = class.names[j]
      temp$length = end(match) - start(match) + 1
      allRepeats = rbind(allRepeats, temp)
    }
  }
  gc()
  time = Sys.time() - time
  print(time)
}

allRepeats$seqSameOrientation = allRepeats$x

for( i in 1 : length( allRepeats$x ))
{
  print(i)
  if(allRepeats$class[i] == "cen180minus")
  {
    allRepeats$seqSameOrientation[i] = revCompString(allRepeats$seqSameOrientation[i])
  }
}

allRepeats$length = allRepeats$end - allRepeats$start + 1

iii = nrow(allRepeats)
while(iii > 1)
{
  overlap = allRepeats$start[iii] - allRepeats$end[iii - 1] - 1
  if(overlap < 0 & overlap >= -10)
  {
    print(overlap)
    overlap = -overlap
    b = overlap %/% 2
    c = overlap - b
    allRepeats$start[iii] = allRepeats$start[iii] + b
    allRepeats$end[iii - 1] = allRepeats$end[iii - 1] - c
    
  } else if(overlap < -10 & overlap > -500)
  {
    if(allRepeats$length[iii] <= allRepeats$length[iii - 1])
    {
      allRepeats = allRepeats[-c(iii),] 
    } else
    {
      allRepeats = allRepeats[-c(iii - 1),] 
    }
  }
  iii = iii - 1
}

write.csv(x = allRepeats, file = "C:/Users/wlodz/Desktop/CEN180allrepeats_secondConsensus_col0601_sameOrientationToAlign.csv", row.names = FALSE)

for(i in 1 : length(allRepeats$x))
{
  #print(i)
  write.fasta(allRepeats$seqSameOrientation[i], names = paste(substr(allRepeats$contig[i],2, 5), allRepeats$class[i], allRepeats$start[i], allRepeats$end[i], sep = "_"), file.out = "C:/Users/wlodz/Desktop/CEN180allrepeats_secondconsensus_col0601_sameOrientationToAlign.fasta", open = "a")
}

###### mafft, consensus, 

#% mafft-sparsecore.rb -A '--kimura 1' -C '--globalpair --maxiterate 3 --kimura 1' -o i -p 1000 -i input 







cen180 = read.csv(file = "C:/Users/wlodz/Desktop/col1227/cen180v0.3.csv")

CEN180alignment = read.alignment(file = "C:/Users/wlodz/Desktop/alignmentsForHORs/col0601.alignment.fasta", format = "fasta", forceToLower = TRUE) 
alignmentVector = vector(mode = "character", length = length(CEN180alignment$seq))

for(i in 1: length(alignmentVector))
{
  print(i)
  alignmentVector[i] = CEN180alignment$seq[[i]][1] %>% strsplit(, split = "")
}
remove(CEN180alignment)
gc()

consensusA = vector(mode = "numeric", length = length(alignmentVector[[1]]))
consensusT = vector(mode = "numeric", length = length(alignmentVector[[1]]))
consensusC = vector(mode = "numeric", length = length(alignmentVector[[1]]))
consensusG = vector(mode = "numeric", length = length(alignmentVector[[1]]))
consensusN = vector(mode = "numeric", length = length(alignmentVector[[1]]))

increment = 1/nrow(cen180)

for(i in 1 : nrow(cen180))
{
  print(i)
  for(j in 1 : length(alignmentVector[[1]]))
  {
    consensusA[j] = consensusA[j] + increment * (alignmentVector[[i]][j] == "a")
    consensusT[j] = consensusT[j] + increment * (alignmentVector[[i]][j] == "t")
    consensusC[j] = consensusC[j] + increment * (alignmentVector[[i]][j] == "c")
    consensusG[j] = consensusG[j] + increment * (alignmentVector[[i]][j] == "g")
    consensusN[j] = consensusN[j] + increment * (alignmentVector[[i]][j] == "-")
  }
}

consensusCEN180 = vector(mode = "character", length = length(alignmentVector[[1]]))
for(i in 1 : length(alignmentVector[[1]]))
{
  top = sort(c(consensusA[i], consensusT[i], consensusC[i], consensusG[i], consensusN[i]), decreasing = TRUE)[1]
  if(top == consensusA[i])
  {
    consensusCEN180[i] = "a"
  } else if(top == consensusC[i])
  {
    consensusCEN180[i] = "c"
  } else if(top == consensusT[i])
  {
    consensusCEN180[i] = "t"
  } else if(top == consensusG[i])
  {
    consensusCEN180[i] = "g"
  } else if(top == consensusN[i])
  {
    consensusCEN180[i] = "-"
  }
}


cen180$weightedSNV = 0

#vecnum = vector(mode = "numeric", length = 748000)
#vecstr = vector(mode = "character", length = 748000)
#SNVcoordinates = data.frame(vecnum,vecnum,vecstr,vecstr)
#names(SNVcoordinates) = c("chromosome", "chromosome_coordinate", "variant_nucleotide", "consensus_nucleotide")
time = Sys.time()
#a = 1
for(i in 1 : nrow(cen180))
{
  print(i)
  #if(cen180$strand[i] == "+")
  #{
    for(j in 1 : length(alignmentVector[[1]]))
    {
      cen180$weightedSNV[i] = cen180$weightedSNV[i] + consensusA[j] * (alignmentVector[[i]][j] != "a") + consensusC[j] * (alignmentVector[[i]][j] != "c") + consensusT[j] * (alignmentVector[[i]][j] != "t") + consensusG[j] * (alignmentVector[[i]][j] != "g") + consensusN[j] * (alignmentVector[[i]][j] != "-") 
      #if(alignmentVector[[i]][j] != consensusCEN180[j])
      #{
      #  SNVcoordinates$chromosome[a] = cen180$chromosome[i]
      #  SNVcoordinates$chromosome_coordinate[a] = cen180$start[i] + length(alignmentVector[[i]][1:j-1][alignmentVector[[i]][1:j-1]!="-"]) 
      #  SNVcoordinates$variant_nucleotide[a] = alignmentVector[[i]][j]
      #  SNVcoordinates$consensus_nucleotide[a] = consensusCEN180[j]
      #  a = a + 1
      #}
    }
  #} else
  #{
   # for(j in 1 : 968)
    #{
     # cen180$weightedSNV[i] = cen180$weightedSNV[i] + consensusA[j] * (alignmentVector[[i]][j] != "a") + consensusC[j] * (alignmentVector[[i]][j] != "c") + consensusT[j] * (alignmentVector[[i]][j] != "t") + consensusG[j] * (alignmentVector[[i]][j] != "g") + consensusN[j] * (alignmentVector[[i]][j] != "-") 
      #if(alignmentVector[[i]][j] != consensusCEN180[j])
      #{
      #  SNVcoordinates$chromosome[a] = cen180$chromosome[i]
      #  SNVcoordinates$chromosome_coordinate[a] = cen180$end[i] - length(alignmentVector[[i]][1:j-1][alignmentVector[[i]][1:j-1]!="-"])
      #  SNVcoordinates$variant_nucleotide[a] = alignmentVector[[i]][j]
      #  SNVcoordinates$consensus_nucleotide[a] = consensusCEN180[j]
      #  a = a + 1
      #}
    #}
  #}
}
print(Sys.time() - time)





# make BED files to import annotations into genious:
# - tab separated values table
# - first column: name of the file
# - second column: start coordinate (counting from 0) (so every coordinate from my table should be made -1)
# - third column: end coordinate
# - fourth column: annotation name
# - fifth column: don't know, just put 0
# - sixth column: strand (+/-)
# afetr exporting, add "track name="No Track"" line at the beginning

firstC = vector(mode = "character", length = length(cen180$x))
secondC = cen180$start
thirdC = cen180$end
fourthC = paste("chr", cen180$chromosome, cen180$start, cen180$end, sep = "_")
fifthC = vector(mode = "numeric", length = length(cen180$x))
sixthC = vector(mode = "numeric", length = length(cen180$x))
for(i in 1 : nrow(cen180))
{
  print(i)
  firstC[i] = "chromosome"
  fifthC[i] = 0
  if(cen180$strand[i] == "-")
  {
    sixthC[i] = "-"
  }  else
  {
    sixthC[i] = "+"
  }
}
bedfile = data.frame(firstC, secondC, thirdC, fourthC, fifthC, sixthC)
bedfile1 = bedfile[substr(bedfile$fourthC,5,5) == "1",]
bedfile2 = bedfile[substr(bedfile$fourthC,5,5) == "2",]
bedfile3 = bedfile[substr(bedfile$fourthC,5,5) == "3",]
bedfile4 = bedfile[substr(bedfile$fourthC,5,5) == "4",]
bedfile5 = bedfile[substr(bedfile$fourthC,5,5) == "5",]

write.table(x = bedfile1, file = "bedannotations1.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(x = bedfile2, file = "bedannotations2.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(x = bedfile3, file = "bedannotations3.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(x = bedfile4, file = "bedannotations4.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(x = bedfile5, file = "bedannotations5.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


cen180.per.chromosome = vector()
for(i in 1 : 5)
{
  cen180.per.chromosome[i] = length(cen180$X[cen180$chromosome == i])
}

write.csv(x = cen180.per.chromosome, file = "0601cen180.per.chromosome.csv")




cen180.genome.scale.unique = 0
cen180.chromosome.scale.unique = vector(length = 5, mode = "numeric")

for(ii in 1 : 5)
{
  print(ii)
  temp = cen180[as.numeric(cen180$chromosome) == ii,]
  for(i in 1 : nrow(temp))
  {
    cen180.chromosome.scale.unique[ii] = cen180.chromosome.scale.unique[ii] + 1 * ((nrow(temp[temp$seq == temp$seq[i],])) > 1)
    #cen180$repeated.chromosome[i] = nrow(temp[temp$seq == temp$seq[i],])
  }
}

for(i in 1 : nrow(cen180))
{
  print(i)
  cen180.genome.scale.unique = cen180.genome.scale.unique + 1 * ((nrow(cen180[cen180$seq == cen180$seq[i],])) > 1)
  #cen180$repeated.genome[i] = nrow(cen180[cen180$seq == cen180$seq[i],])
}


cen180.stats = data.frame(average.SNV.to.consensus = ave(cen180$weightedSNV)[1], average.length = ave(cen180$length)[1], average.HOR.count = ave(cen180$HORcount)[1], 
                          average.HOR.lengths.sum = ave(cen180$HORlengthsSum)[1], percentage.repeats.with.no.HORs = (100 * nrow(cen180[cen180$HORcount == 0,])/ nrow(cen180)))

cen180$repeated.chromosome = cen180$repeated.chromosome - 1
  
cen180$repeated.genome = cen180$repeated.genome - 1
  
write.csv(x = cen180.chromosome.scale.unique, file = "0601cen180.chromosome.scale.unique.csv")
write.csv(x = cen180.stats, file = "0601cen180.stats.csv")


write.csv(x = cen180, file = "0601cen180V0.6.csv")







cen180$chromosome = as.numeric(str_sub(cen180$contig, 5, 5))


hors = read.csv(file= "C:/Users/wlodz/Desktop/alignmentsForHORs/col0601.alignment.fasta_t_5_c_2.csv", header = TRUE)
hors = as.data.frame(hors)



#cen180 = read.csv("C:/Users/wlodz/Desktop/col1227/cen180v0.3.csv", header = TRUE)
#cen180 = cen180[,-c(1,4,9,11,12)]

hors$start.A.bp = cen180$start[hors$START.A]
hors$start.B.bp = cen180$start[hors$START.B]
hors$end.A.bp = cen180$end[hors$END.A]
hors$end.B.bp = cen180$end[hors$END.B] 
hors$chrA = cen180$chromosome[hors$START.A]  
hors$chrB = cen180$chromosome[hors$START.B]
#hors$adjusted.start.A.bp =  cen180$adjustedStart[hors$START.A]
#hors$adjusted.start.B.bp =  cen180$adjustedStart[hors$START.B]  
hors$block.size.in.units = hors$END.A - hors$START.A + 1
hors$block.A.size.bp = hors$end.A.bp - hors$start.A.bp + 1
hors$block.B.size.bp = hors$end.B.bp - hors$start.B.bp + 1
hors = hors[hors$block.size.in.units > 2]

##### check if HOR is spanning a gap and remove those that do
#horsB = hors
for(i in 1 : length(hors$index))
{
  block1size =abs(hors$end.A.bp[i] - hors$start.A.bp[i])
  block2size = abs(hors$end.B.bp[i] - hors$start.B.bp[i])
  if(abs(block2size - block1size) > 1000)
  {
    print(i)
    hors = hors[-i,]
  }
}
#removed 26 HORs
#remove(horsB)
##############################################################


# check overlapping HORs and remove them
elo = 0
for( i in 1 : length(hors$index))
{
  if(hors$DIR.2.same_1.opposite.[i] == 2)
  {
    if(hors$START.B[i] - hors$END.A[i] < 1)
    {
      elo = c(elo, i)
      #hors = hors[-i,]
    }
  }
}
elo = elo[-1]
#2332 cases
horsb = hors
hors = hors[-elo,]
remove(horsb)

elo = 0
for( i in 1 : length(hors$index))
{
  if(hors$DIR.2.same_1.opposite.[i] == 1)
  {
    if(hors$START.B[i] - hors$END.A[i] < 1)
    {
      elo = c(elo, i)
    }
  }
}
#0 cases

hors = hors[(hors$block.A.size.bp - (hors$block.size.in.units * 178) <= 100) & (hors$block.B.size.bp - (hors$block.size.in.units * 178) <= 100),]


length(hors$index[hors$DIR.2.same_1.opposite.==2])
#2204518
length(hors$index[hors$DIR.2.same_1.opposite.==1])
#1624178
###########################
cen180backup = cen180

cen180$start[cen180$chromosome == "1"] = cen180$start[cen180$chromosome == "1"] + 14328593
cen180$end[cen180$chromosome == "1"] = cen180$end[cen180$chromosome == "1"] + 14328593

cen180$start[cen180$chromosome == "2"] = cen180$start[cen180$chromosome == "2"] + 3334745
cen180$end[cen180$chromosome == "2"] = cen180$end[cen180$chromosome == "2"] + 3334745

cen180$start[cen180$chromosome == "3"] = cen180$start[cen180$chromosome == "3"] + 12820378
cen180$end[cen180$chromosome == "3"] = cen180$end[cen180$chromosome == "3"] + 12820378

cen180$start[cen180$chromosome == "4"] = cen180$start[cen180$chromosome == "4"] + 2988897
cen180$end[cen180$chromosome == "4"] = cen180$end[cen180$chromosome == "4"] + 2988897

cen180$start[cen180$chromosome == "5"] = cen180$start[cen180$chromosome == "5"] + 10996544
cen180$end[cen180$chromosome == "5"] = cen180$end[cen180$chromosome == "5"] + 10996544

##############################

#plot nicely

cen.rough.start <- c(14590750, 3474530, 13347090, 2773972, 11783865 - 250000)
cen.rough.end <- c(17808182, 6196091, 15983029, 7226978, 14551874 + 250000)

png(filename = paste("col0601MNew HOR 15 plots", ".png", sep = ""), width = 10000, height = 10000, pointsize = 50)
par(mfrow = c(5,5))
for(j in 5 : 1)
{
  for(i in 1 : 5)
  {
    if(i > j)
    {
      plot.new()
    } else
    {
      print(i)
      plot(xlab = paste("chromosome ", i, ", bp", sep = ""), ylab = paste("chromosome ", j, ", bp", sep = ""), x = hors$start.A.bp[(hors$chrA == i & hors$chrB == j) | (hors$chrA == j & hors$chrB == i)], pch = 19, cex = 0.1, y = hors$start.B.bp[(hors$chrA == i & hors$chrB == j) | (hors$chrA == j & hors$chrB == i)], xlim = c(cen.rough.start[i], cen.rough.end[i]), ylim = c(cen.rough.start[j], cen.rough.end[j]))
    }
  }
}
dev.off()




#count HORs


# make activity (count x size) plots

cen180$HORcount = 0
cen180$HORlengthsSum = 0

for(i in 1 : length(hors$index))
{
  print(i)
  for(j in 0 : (hors$block.size.in.units[i] - 1))
  {
    cen180$HORcount[hors$START.A[i] + j] = cen180$HORcount[hors$START.A[i] + j] + 1
    cen180$HORcount[hors$START.B[i] + j] = cen180$HORcount[hors$START.B[i] + j] + 1
    cen180$HORlengthsSum[hors$START.A[i] + j] = cen180$HORlengthsSum[hors$START.A[i] + j] + hors$block.size.in.units[i]
    cen180$HORlengthsSum[hors$START.B[i] + j] = cen180$HORlengthsSum[hors$START.B[i] + j] + hors$block.size.in.units[i]
  }
}


cen180$averagedSUMofHORlengths = 0
for(i in 25 : (nrow(cen180) - 25))
{
  cen180$averagedSUMofHORlengths[i] = ave(c(cen180$HORlengthsSum[i - 24]:cen180$HORlengthsSum[i + 24]))[1]
}

png(filename = "0601MCEN180 activity new.png", width = 2000, height = 2500, pointsize = 25)
par(mfrow=c(5,1))
plot(cen180$start[cen180$chromosome == 1], cen180$averagedSUMofHORlengths[cen180$chromosome == 1], cex.lab = 1.5, cex.axis = 1.5, cex = 0.5, pch = 19, ylim = c(0,2000), main ="Chromosome 1 repeats. Sum of HORs sizes each repeat is a part of", ylab = "Activity score", xlab = "Position, bp")
plot(cen180$start[cen180$chromosome == 2], cen180$averagedSUMofHORlengths[cen180$chromosome == 2], cex.lab = 1.5, cex.axis = 1.5, cex = 0.5, pch = 19, ylim = c(0,2000),  main ="Chromosome 2 repeats. Sum of HORs sizes each repeat is a part of", ylab = "Activity score", xlab = "Position, bp")
plot(cen180$start[cen180$chromosome == 3], cen180$averagedSUMofHORlengths[cen180$chromosome == 3], cex.lab = 1.5, cex.axis = 1.5, cex = 0.5, pch = 19, ylim = c(0,2000),  main ="Chromosome 3 repeats. Sum of HORs sizes each repeat is a part of", ylab = "Activity score", xlab = "Position, bp")
plot(cen180$start[cen180$chromosome == 4], cen180$averagedSUMofHORlengths[cen180$chromosome == 4], cex.lab = 1.5, cex.axis = 1.5, cex = 0.5, pch = 19, ylim = c(0,2000),  main ="Chromosome 4 repeats. Sum of HORs sizes each repeat is a part of", ylab = "Activity score", xlab = "Position, bp")
plot(cen180$start[cen180$chromosome == 5], cen180$averagedSUMofHORlengths[cen180$chromosome == 5], cex.lab = 1.5, cex.axis = 1.5, cex = 0.5, pch = 19, ylim = c(0,2000),  main ="Chromosome 5 repeats. Sum of HORs sizes each repeat is a part of", ylab = "Activity score", xlab = "Position, bp")
dev.off()


#### export data

write.csv(x = hors, file = "HORsV0.4.csv", row.names = FALSE)

write.csv(x = cen180, file = "cen180sV0.5.csv", row.names = FALSE)
#hors.nomonomer = hors[,-c(1,2,3,4,5,14,15,16)]
#write.csv(x = hors.nomonomer, file = "HORsV0.4.no.monomer.info.csv", row.names = FALSE)


#cen.rough.start <- c(14590750, 3474530, 13347090, 2773972, 11783865 - 250000)
#cen.rough.end <- c(17808182, 6196091, 15983029, 7226978, 14551874 + 250000)


#for(i in 1 : 20)
#{
#  print(i)
#  for(j in 1 : 5)
#  {
#    print(j)
#    for(k in 1 : 2)
#    {
#      horsChr = hors[hors$chrA == j & hors$chrB == j & hors$SNVcount/hors$block.size.in.units <= i/2 & hors$DIR.2.same_1.opposite. == k,]
#      png(filename = paste("HORplots/chr ", j, " HORs with ", k, " strand heatmap.png", sep = ""), width = 1600, height = 1500, pointsize = 35)
#      print(ggplot(data = data.frame(startA = horsChr$start.A.bp, startB = horsChr$start.B.bp), aes(x = startA, y = startB)) + geom_bin2d(bins = 1000) + scale_fill_gradientn(limits=c(0,70), breaks=seq(0, 60, by=10), colours=rainbow(6)) + ggtitle(paste("chr ", j, " HORs with ", k, " strand heatmap", sep = "")) + theme(plot.title = element_text(size = 50, face="bold")) + xlim(cen.rough.start[j], cen.rough.end[j]) + ylim(cen.rough.start[j], cen.rough.end[j]))
#      dev.off()
#    }
#  }
#}





#j = 5
#length(horse$index[horse$chrA == j & horse$chrB == j])
#length(horse$index[horse$chrA == j & horse$chrB == j & horse$DIR.2.same_1.opposite. == 2])
#length(horse$index[horse$chrA == j & horse$chrB == j & horse$DIR.2.same_1.opposite. == 1])
#
#i = 2
#for(j in 1 : 5)
#{
#  print(j)
#  for(k in 1 : 2)
#  {
#    horsChr = hors[hors$chrA == j & hors$chrB == j & hors$block.size.in.units >= i & hors$DIR.2.same_1.opposite. == k,]
#    png(filename = paste("HORplots/Ochr ", j, " HORs with ", k, " strand heatmap.png", sep = ""), width = 1600, height = 1500, pointsize = 35)
#    print(ggplot(data = data.frame(startA = horsChr$start.A.bp, startB = horsChr$start.B.bp), aes(x = startA, y = startB)) + geom_bin2d(bins = 1000) + scale_fill_gradientn(limits=c(0,70), breaks=seq(0, 60, by=10), colours=rainbow(6)) + ggtitle(paste("chr ", j, " HORs with ", k, " strand heatmap", sep = "")) + theme(plot.title = element_text(size = 50, face="bold")) + xlim(cen.rough.start[j], cen.rough.end[j]) + ylim(cen.rough.start[j], cen.rough.end[j]))
#    dev.off()
#  }
#}

hor.matrix = matrix(data = 0 ,nrow = 5, ncol = 5, dimnames = list(c("chr1", "chr2", "chr3", "chr4", "chr5"),c("chr1", "chr2", "chr3", "chr4", "chr5")))
for(i in 1 : 5)
{
  for(ii in i : 5)
  {
    hor.matrix[i,ii] = length(hors$index[hors$chrA == i & hors$chrB == ii])
  }
}

write.csv(x = hor.matrix, file = "0601hor.matrix.csv")


hors.stats = data.frame(average.block.size.in.units = ave(hors$block.size.in.units)[1], biggest.block.size.in.units = max(hors$block.size.in.units), 
                        average.hor.blocks.distance.in.bp = ave(abs(hors$start.B.bp[hors$chrA == hors$chrB] - hors$start.A.bp[hors$chrA == hors$chrB]))[1],
                        max.hor.blocks.distance.in.bp = max(abs(hors$start.B.bp[hors$chrA == hors$chrB] - hors$start.A.bp[hors$chrA == hors$chrB])))

hors[which.max(hors$block.size.in.units),]



write.csv(x = hors.stats, file = "0601hors.stats.csv")



####CEN160






cen160plus = tolower("GTCAAATGCATTGGATTGTGACACATTTTGACCATAGAAACACTAACAAAGCTATTTACTGCTTCTAAGCAATTTTTTGTTGGTTTTAGCCTCTTTTGGGAGAAAATGGGTATAAGTGTTGTCTAAACACTCCTAATCCATCTCTAACTCTTATAATTA")
cen160minus = tolower("TAATTATAAGAGTTAGAGATGGATTAGGAGTGTTTAGACAACACTTATACCCATTTTCTCCCAAAAGAGGCTAAAACCAACAAAAAATTGCTTAGAAGCAGTAAATAGCTTTGTTAGTGTTTCTATGGTCAAAATGTGTCACAATCCAATGCATTTGAC")

class.names = c("cen160plus", "cen160minus")

class.sequences = c(cen160plus, cen160minus)

allRepeats = data.frame()
allRepeats$contig = vector(mode = "numeric")
allRepeats$class = vector(mode = "numeric")
allRepeats$start = vector(mode = "numeric")
allRepeats$end = vector(mode = "numeric")
allRepeats$sequence = vector(mode = "character")
allRepeats$strand = vector(mode = "character")

a = read.fasta(file = "C:/Users/wlodz/Desktop/col1227/t2t-col.20201227.fasta", as.string = TRUE)


for(i in 1 : 1)
{
  time = Sys.time()
  print(paste("fasta file ", i, "/7", sep = ""))
  for(j in 1 : 2)#length(class.names))
  {
    print(paste("repeat class: ", class.names[j], sep = ""))
    maxMis = nchar(class.sequences[1]) %/% 3
    match = matchPattern(pattern = class.sequences[j], subject = a[[i]][1], max.mismatch = maxMis, with.indels = TRUE)
    temp = as.data.frame(match)
    if(length(temp$x) > 0)
    {
      print(i)
      print(j)
      temp$contig = getAnnot(a[[i]])
      temp$start  = start(match)
      temp$end = end(match)
      temp$strand = '+'
      temp$class = class.names[j]
      temp$length = end(match) - start(match) + 1
      allRepeats = rbind(allRepeats, temp)
    }
  }
  gc()
  time = Sys.time() - time
  print(time)
}

allRepeats$seqSameOrientation = allRepeats$x

for( i in 1 : length( allRepeats$x ))
{
  print(i)
  if(allRepeats$class[i] == "cen160minus")
  {
    allRepeats$seqSameOrientation[i] = revCompString(allRepeats$seqSameOrientation[i])
  }
}

allRepeats$length = allRepeats$end - allRepeats$start + 1

iii = nrow(allRepeats)
while(iii > 1)
{
  overlap = allRepeats$start[iii] - allRepeats$end[iii - 1] - 1
  if(overlap < 0 & overlap >= -10)
  {
    print(overlap)
    overlap = -overlap
    b = overlap %/% 2
    c = overlap - b
    allRepeats$start[iii] = allRepeats$start[iii] + b
    allRepeats$end[iii - 1] = allRepeats$end[iii - 1] - c
    
  } else if(overlap < -10 & overlap > -500)
  {
    if(allRepeats$length[iii] <= allRepeats$length[iii - 1])
    {
      allRepeats = allRepeats[-c(iii),] 
    } else
    {
      allRepeats = allRepeats[-c(iii - 1),] 
    }
  }
  iii = iii - 1
}

write.csv(x = allRepeats, file = "C:/Users/wlodz/Desktop/col1227/CEN160allrepeat.csv", row.names = FALSE)



#add GC to frame6
frame6 = read.csv(file = "C:/Users/wlodz/Desktop/mafft-7.475-win64-signed/frame6.csv", header = TRUE)


for(i in 1 : nrow(frame6))
{
  split = strsplit(frame6$consensus.sequence[i], split = "")[[1]]
  frame6$GC[i] = GC(split)
  frame6$G[i] = length(which(split == "G")) / length(split)
  frame6$C[i] = length(which(split == "C")) / length(split)
  frame6$A[i] = length(which(split == "A")) / length(split)
  frame6$T[i] = length(which(split == "T")) / length(split)
}
