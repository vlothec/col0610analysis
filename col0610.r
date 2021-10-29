library(Biostrings)
library(stringr)
library(beepr)
library(svMisc)
library(seqinr)
library(base)
library(rstudioapi)
library(magrittr)
library(grDevices)
revCompString = function(DNAstr) {  return(tolower(toString(reverseComplement(DNAString(DNAstr))))) }
library(ggplot2)
library(ggseqlogo)



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

a = read.fasta(file = "C:/Users/wlodz/Desktop/t2t-col.20210610.fasta", as.string = TRUE)


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

chr1 = allRepeats[allRepeats$contig == ">Chr1",]
chr1 = chr1[order(chr1$start),]

chr2 = allRepeats[allRepeats$contig == ">Chr2",]
chr2 = chr2[order(chr2$start),]

chr3 = allRepeats[allRepeats$contig == ">Chr3",]
chr3 = chr3[order(chr3$start),]

chr4 = allRepeats[allRepeats$contig == ">Chr4",]
chr4 = chr4[order(chr4$start),]

chr5 = allRepeats[allRepeats$contig == ">Chr5",]
chr5 = chr5[order(chr5$start),]

#allRepeatsb = allRepeats
allRepeats = rbind(chr1,chr2,chr3,chr4,chr5)

write.csv(x = allRepeats, file = "C:/Users/wlodz/Desktop/CEN180allrepeats_T2T0610_sameOrientationToAlign.csv", row.names = FALSE)

for(i in 1 : length(allRepeats$x))
{
  #print(i)
  write.fasta(allRepeats$seqSameOrientation[i], names = paste(substr(allRepeats$contig[i],2, 5), allRepeats$class[i], allRepeats$start[i], allRepeats$end[i], sep = "_"), file.out = "C:/Users/wlodz/Desktop/CEN180allrepeats_T2T0610_sameOrientationToAlign.fasta", open = "a")
}

###### mafft, consensus, 

#% mafft-sparsecore.rb -A '--kimura 1' -C '--globalpair --maxiterate 3 --kimura 1' -o i -p 1000 -i input 

#mafft --thread 8 --threadtb 5 --threadit 0 --inputorder --kimura 1 --large --retree 2 input > output





cen180 = allRepeats

CEN180alignment = read.alignment(file = "C:/Users/wlodz/Desktop/alignmentsForHORs/T2T0610.alignment.fasta", format = "fasta", forceToLower = TRUE) 

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

increment = 1/length(alignmentVector)

for(i in 1 : length(alignmentVector))
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

time = Sys.time()

for(i in 1 : nrow(cen180))
{
  print(i)
  
  for(j in 1 : length(alignmentVector[[1]]))
  {
    cen180$weightedSNV[i] = cen180$weightedSNV[i] + consensusA[j] * (alignmentVector[[i]][j] != "a") + consensusC[j] * (alignmentVector[[i]][j] != "c") + consensusT[j] * (alignmentVector[[i]][j] != "t") + consensusG[j] * (alignmentVector[[i]][j] != "g") + consensusN[j] * (alignmentVector[[i]][j] != "-") 
    
  }
  
}
print(Sys.time() - time)



consensusT2T061010 = data.frame(consensusA, consensusT, consensusC, consensusG, consensusN)
write.csv(consensusT2T061010, file = "C:/Users/wlodz/Desktop/consensusT2T061010.csv")


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
  firstC[i] = paste("chromosome", cen180$chromosome[i])
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

write.table(x = bedfile1, file = "C:/Users/wlodz/Desktop/bedannotations1T2T0610.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(x = bedfile2, file = "C:/Users/wlodz/Desktop/bedannotations2T2T0610.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(x = bedfile3, file = "C:/Users/wlodz/Desktop/bedannotations3T2T0610.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(x = bedfile4, file = "C:/Users/wlodz/Desktop/bedannotations4T2T0610.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(x = bedfile5, file = "C:/Users/wlodz/Desktop/bedannotations5T2T0610.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


cen180.per.chromosome = vector()
for(i in 1 : 5)
{
  cen180.per.chromosome[i] = length(cen180$X[cen180$chromosome == i])
}

write.csv(x = cen180.per.chromosome, file = "C:/Users/wlodz/Desktop/T2T0610.per.chromosome.csv")

write.csv(x = cen180, file = "C:/Users/wlodz/Desktop/T2T0610.v.1.0.csv")

###########################
#                         #
#  MAFFT alignment here   #
#                         #
##################### DO HORs T 5 C 1


cen180$chromosome = as.numeric(str_sub(cen180$contig, 5, 5)) #TODO


hors = read.csv(file= "C:/Users/wlodz/Desktop/alignmentsForHORs/T2T0610.alignment.fasta_t_5_c_1.csv", header = TRUE)
hors = as.data.frame(hors)



hors$start.A.bp = cen180$start[hors$START.A]
hors$start.B.bp = cen180$start[hors$START.B]
hors$end.A.bp = cen180$end[hors$END.A]
hors$end.B.bp = cen180$end[hors$END.B] 
hors$chrA = cen180$chromosome[hors$START.A]  
hors$chrB = cen180$chromosome[hors$START.B]
hors$block.size.in.units = hors$END.A - hors$START.A + 1
hors$block.A.size.bp = hors$end.A.bp - hors$start.A.bp + 1
hors$block.B.size.bp = hors$end.B.bp - hors$start.B.bp + 1
#hors = hors[hors$block.size.in.units > 2,]

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
#1568 cases
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


#plot nicely

cen.rough.start <- c(14590750, 3474530, 13347090, 2773972, 11783865 - 250000)
cen.rough.end <- c(17808182, 6196091, 15983029, 7226978, 14551874 + 250000)

png(filename = paste("C:/Users/wlodz/Desktop/T2T0610 HOR 15 plots", ".png", sep = ""), width = 10000, height = 10000, pointsize = 50)
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

png(filename = "C:/Users/wlodz/Desktop/T2T0610 activity new.png", width = 2000, height = 2500, pointsize = 25)
par(mfrow=c(5,1))
plot(cen180$start[cen180$chromosome == 1], cen180$averagedSUMofHORlengths[cen180$chromosome == 1], cex.lab = 1.5, cex.axis = 1.5, cex = 0.5, pch = 19, ylim = c(0,max(cen180$averagedSUMofHORlengths)+100), main ="Chromosome 1 repeats. Sum of HORs sizes each repeat is a part of", ylab = "Activity score", xlab = "Position, bp")
plot(cen180$start[cen180$chromosome == 2], cen180$averagedSUMofHORlengths[cen180$chromosome == 2], cex.lab = 1.5, cex.axis = 1.5, cex = 0.5, pch = 19, ylim = c(0,max(cen180$averagedSUMofHORlengths)+100),  main ="Chromosome 2 repeats. Sum of HORs sizes each repeat is a part of", ylab = "Activity score", xlab = "Position, bp")
plot(cen180$start[cen180$chromosome == 3], cen180$averagedSUMofHORlengths[cen180$chromosome == 3], cex.lab = 1.5, cex.axis = 1.5, cex = 0.5, pch = 19, ylim = c(0,max(cen180$averagedSUMofHORlengths)+100),  main ="Chromosome 3 repeats. Sum of HORs sizes each repeat is a part of", ylab = "Activity score", xlab = "Position, bp")
plot(cen180$start[cen180$chromosome == 4], cen180$averagedSUMofHORlengths[cen180$chromosome == 4], cex.lab = 1.5, cex.axis = 1.5, cex = 0.5, pch = 19, ylim = c(0,max(cen180$averagedSUMofHORlengths)+100),  main ="Chromosome 4 repeats. Sum of HORs sizes each repeat is a part of", ylab = "Activity score", xlab = "Position, bp")
plot(cen180$start[cen180$chromosome == 5], cen180$averagedSUMofHORlengths[cen180$chromosome == 5], cex.lab = 1.5, cex.axis = 1.5, cex = 0.5, pch = 19, ylim = c(0,max(cen180$averagedSUMofHORlengths)+100),  main ="Chromosome 5 repeats. Sum of HORs sizes each repeat is a part of", ylab = "Activity score", xlab = "Position, bp")
dev.off()


#### export data

write.csv(x = hors, file = "T2T0610HORsV1.0.csv", row.names = FALSE)

cen180$percentageIdentity = 100-(cen180$weightedSNV*100/809)#TODO
write.csv(x = cen180, file = "T2T0610cen180sV.1.0.csv", row.names = FALSE)



hor.matrix = matrix(data = 0 ,nrow = 5, ncol = 5, dimnames = list(c("chr1", "chr2", "chr3", "chr4", "chr5"),c("chr1", "chr2", "chr3", "chr4", "chr5")))
for(i in 1 : 5)
{
  for(ii in i : 5)
  {
    hor.matrix[i,ii] = length(hors$index[hors$chrA == i & hors$chrB == ii])
  }
}

write.csv(x = hor.matrix, file = "T2T0610v.v.0.hor.matrix.csv")

##########################################################

cen180 = read.csv("C:/Users/wlodz/Desktop/T2T0610cen180sV.1.0.csv")

consensusCEN180 = tolower(strsplit("--------------------------AGTA-T-AA-G-AA-CTT---A-A-A-CC-GC------AA-----------C-C-------C-G--A------T-C-T---------T----AAA--------A-G--------------CC-T-------AAG------TA-G--------T------GT-TT-------------C-CTTGT----T-------A-G----------------------------A-------AGAC-AC---------AA------AG-C-----------C----AAA----------G----------A---------------------------------CT-------------------C-ATA---------------------------------------T-G-GAC-TTT----------------------------GGCTAC-A-CC---------A--T-----G--AAAG-----------C-----TT-T-G--------------------AG------AA----GC-----------------AA-G-A--------------------A--G----AA-G------------G-T-TG--------------------G-TT-A------G-------TG-T-T-TT-G-GA-----------G-TCG-----A----A------TA-T-------------GAC-------T-TG--------A-----TG-T-----CAT------GTGTA-TGATTG------------", split = "")[[1]])
vecnum = vector(mode = "numeric", length = 1000000)
vecstr = vector(mode = "character", length = 1000000)
SNVcoordinates = data.frame(vecnum,vecnum,vecstr,vecstr)
names(SNVcoordinates) = c("chromosome", "chromosome_coordinate", "variant_nucleotide", "consensus_nucleotide")
time = Sys.time()
a = 1
for(i in 1 : length(alignmentVector))
{
  print(i)
  if(cen180$strand.1[i] == "cen180plus")
  {
    for(j in 1 : length(alignmentVector[[1]]))
    {
      #cen180$weightedSNV[i] = cen180$weightedSNV[i] + consensusA[j] * (alignmentVector[[i]][j] != "a") + consensusC[j] * (alignmentVector[[i]][j] != "c") + consensusT[j] * (alignmentVector[[i]][j] != "t") + consensusG[j] * (alignmentVector[[i]][j] != "g") + consensusN[j] * (alignmentVector[[i]][j] != "-") 
      if(alignmentVector[[i]][j] != consensusCEN180[j])
      {
        SNVcoordinates$chromosome[a] = cen180$chromosome[i]
        SNVcoordinates$chromosome_coordinate[a] = cen180$start[i] + length(alignmentVector[[i]][1:j-1][alignmentVector[[i]][1:j-1]!="-"]) 
        SNVcoordinates$variant_nucleotide[a] = alignmentVector[[i]][j]
        SNVcoordinates$consensus_nucleotide[a] = consensusCEN180[j]
        a = a + 1
      }
    }
  } else
  {
    for(j in 1 : length(alignmentVector[[1]]))
    {
      #cen180$weightedSNV[i] = cen180$weightedSNV[i] + consensusA[j] * (alignmentVector[[i]][j] != "a") + consensusC[j] * (alignmentVector[[i]][j] != "c") + consensusT[j] * (alignmentVector[[i]][j] != "t") + consensusG[j] * (alignmentVector[[i]][j] != "g") + consensusN[j] * (alignmentVector[[i]][j] != "-") 
      if(alignmentVector[[i]][j] != consensusCEN180[j])
      {
        SNVcoordinates$chromosome[a] = cen180$chromosome[i]
        SNVcoordinates$chromosome_coordinate[a] = cen180$end[i] - length(alignmentVector[[i]][1:j-1][alignmentVector[[i]][1:j-1]!="-"])
        SNVcoordinates$variant_nucleotide[a] = alignmentVector[[i]][j]
        SNVcoordinates$consensus_nucleotide[a] = consensusCEN180[j]
        a = a + 1
      }
    }
  }
}
print(Sys.time() - time)

bSNVcoordinates = SNVcoordinates
SNVcoordinates = SNVcoordinates[1:which(SNVcoordinates$variant_nucleotide == "")[[1]] - 1,]

write.csv(SNVcoordinates, "T2Tcol0610.SNVcoordinates.csv")






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

a = read.fasta(file = "C:/Users/wlodz/Desktop/t2t-col.20210610.fasta", as.string = TRUE)


for(i in 1 : 2)
{
  time = Sys.time()
  print(paste("fasta file ", i, "/7", sep = ""))
  for(j in 1 : length(class.names))
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

write.csv(x = allRepeats, file = "C:/Users/wlodz/Desktop/col0610CEN160allrepeat.csv", row.names = FALSE)


##############chromsome consensus


CEN180alignment = read.alignment(file = "C:/Users/wlodz/Desktop/alignmentsForHORs/T2T0610.alignment.fasta", format = "fasta", forceToLower = TRUE) 

alignmentVector = vector(mode = "character", length = length(CEN180alignment$seq))

for(i in 1: length(alignmentVector))
{
  print(i)
  alignmentVector[i] = CEN180alignment$seq[[i]][1] %>% strsplit(, split = "")
}
remove(CEN180alignment)
gc()

#alignmentVectorb = alignmentVector
#alignmentVector = alignmentVectorb[which(cen180$chromosome == 1)]
#alignmentVector = alignmentVectorb[which(cen180$chromosome == 2)]
#alignmentVector = alignmentVectorb[which(cen180$chromosome == 3)]
#alignmentVector = alignmentVectorb[which(cen180$chromosome == 4)]
#alignmentVector = alignmentVectorb[which(cen180$chromosome == 5)]


consensusA = vector(mode = "numeric", length = length(alignmentVector[[1]]))
consensusT = vector(mode = "numeric", length = length(alignmentVector[[1]]))
consensusC = vector(mode = "numeric", length = length(alignmentVector[[1]]))
consensusG = vector(mode = "numeric", length = length(alignmentVector[[1]]))
consensusN = vector(mode = "numeric", length = length(alignmentVector[[1]]))

increment = 1/length(alignmentVector)

for(i in 1 : length(alignmentVector))
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


cen180.frequencies.chr = data.frame(consensusA, consensusC, consensusG, consensusT)



######### seqLogo
cen180 = read.csv("C:/Users/wlodz/Desktop/T2T0610cen180sV.1.0.csv")
CEN180alignment = read.alignment(file = "C:/Users/wlodz/Desktop/alignmentsForHORs/T2T0610.alignment.fasta", format = "fasta", forceToLower = TRUE) 

alignmentVector = vector(mode = "character", length = length(CEN180alignment$seq))

for(i in 1: length(alignmentVector))
{
  print(i)
  alignmentVector[i] = CEN180alignment$seq[[i]][1] %>% strsplit(, split = "")
}
remove(CEN180alignment)
gc()

alignmentVectorb = alignmentVector

#alignmentVector = alignmentVectorb[which(cen180$chromosome == 1)]
#alignmentVector = alignmentVectorb[which(cen180$chromosome == 2)]
#alignmentVector = alignmentVectorb[which(cen180$chromosome == 3)]
#alignmentVector = alignmentVectorb[which(cen180$chromosome == 4)]
alignmentVector = alignmentVectorb[which(cen180$chromosome == 5)]


consensusA = vector(mode = "numeric", length = length(alignmentVector[[1]]))
consensusT = vector(mode = "numeric", length = length(alignmentVector[[1]]))
consensusC = vector(mode = "numeric", length = length(alignmentVector[[1]]))
consensusG = vector(mode = "numeric", length = length(alignmentVector[[1]]))
consensusN = vector(mode = "numeric", length = length(alignmentVector[[1]]))

increment = 1

for(i in 1 : length(alignmentVector))
{
  for(j in 1 : length(alignmentVector[[1]]))
  {
    consensusA[j] = consensusA[j] + increment * (alignmentVector[[i]][j] == "a")
    consensusT[j] = consensusT[j] + increment * (alignmentVector[[i]][j] == "t")
    consensusC[j] = consensusC[j] + increment * (alignmentVector[[i]][j] == "c")
    consensusG[j] = consensusG[j] + increment * (alignmentVector[[i]][j] == "g")
    consensusN[j] = consensusN[j] + increment * (alignmentVector[[i]][j] == "-")
  }
}

cen180.frequencies.chr = data.frame(consensusA, consensusC, consensusG, consensusT, consensusN)
cen180.frequencies.chr = cen180.frequencies.chr[cen180.frequencies.chr$consensusN < length(alignmentVector)/2,]
cen180.flipped = data.frame(t(cen180.frequencies.chr))
row.names(cen180.flipped) = c("A", "C", "G", "T", "N")


pdf(file = "seqLogo.cen180.chr5.pdf", width = 50, height = 3, pointsize = 20)
#png(filename = "seqLogo.cen180.chr5.png", width = 4000, height = 300)
ggplot() + geom_logo(as.matrix(cen180.flipped), col_scheme='nucleotide') + theme_logo() + theme(axis.text.x = element_text(angle = 90, size=20), axis.text.y = element_text(size=40), axis.title.y = element_text(size=40))
dev.off()




#########################

#cen180 stats unique


cen180$occuredGenome = 0
cen180$occuredChr1 = 0
cen180$occuredChr2 = 0
cen180$occuredChr3 = 0
cen180$occuredChr4 = 0
cen180$occuredChr5 = 0

for(i in 1 : nrow(cen180))
{
  print(i)
  cen180$occuredGenome[i] = length(which(cen180$seqSameOrientation == cen180$seqSameOrientation[i])) - 1
  cen180$occuredChr1[i] = length(which(cen180$seqSameOrientation[cen180$chromosome == 1] == cen180$seqSameOrientation[i]))
  cen180$occuredChr2[i] = length(which(cen180$seqSameOrientation[cen180$chromosome == 2] == cen180$seqSameOrientation[i]))
  cen180$occuredChr3[i] = length(which(cen180$seqSameOrientation[cen180$chromosome == 3] == cen180$seqSameOrientation[i]))
  cen180$occuredChr4[i] = length(which(cen180$seqSameOrientation[cen180$chromosome == 4] == cen180$seqSameOrientation[i]))
  cen180$occuredChr5[i] = length(which(cen180$seqSameOrientation[cen180$chromosome == 5] == cen180$seqSameOrientation[i]))
  if(cen180$chromosome[i] == 1)
  {
    cen180$occuredChr1[i] = cen180$occuredChr1[i] - 1
  } else if(cen180$chromosome[i] == 2)
  {
    cen180$occuredChr2[i] = cen180$occuredChr2[i] - 1
  } else if(cen180$chromosome[i] == 3)
  {
    cen180$occuredChr3[i] = cen180$occuredChr3[i] - 1
  } else if(cen180$chromosome[i] == 4)
  {
    cen180$occuredChr4[i] = cen180$occuredChr4[i] - 1
  } else if(cen180$chromosome[i] == 5)
  {
    cen180$occuredChr5[i] = cen180$occuredChr5[i] - 1
  }
}


cen180.non.chromosome.unique = 0
for(i in 1 : nrow(cen180))
{
  if(cen180[i, 14] + cen180[i, 15] + cen180[i, 16] + cen180[i, 17] + cen180[i, 18] > cen180[i, 13 + cen180$chromosome[i]])
  {
    cen180.non.chromosome.unique = cen180.non.chromosome.unique + 1
  }
}



length(which(duplicated(cen180$seqSameOrientation) | duplicated(cen180$seqSameOrientation, fromLast = TRUE)))

#length(unique(cen180$seqSameOrientation))
length(which(cen180$occuredGenome > 0))

#ratios of internally repeated cen180 sequences in each chromosome
length(which(cen180.1$occuredChr1 > 0)) / nrow(cen180.1)
length(which(cen180.2$occuredChr2 > 0)) / nrow(cen180.2)
length(which(cen180.3$occuredChr3 > 0)) / nrow(cen180.3)
length(which(cen180.4$occuredChr4 > 0)) / nrow(cen180.4)
length(which(cen180.5$occuredChr5 > 0)) / nrow(cen180.5)


#check HORs for the manuscript

HORs = read.csv(file = "C:/Users/wlodz/Desktop/P/T2Tcol0610V1.0.data/T2T0610HORsV1.0.csv")

HORs.same = HORs[HORs$chrA == HORs$chrB,]

HORs.same$distance = abs(HORs.same$start.B.bp - HORs.same$start.A.bp)

nrow(HORs.same[HORs.same$distance < 100000,])/nrow(HORs)

ave(HORs.same$distance)[1]

max(HORs.same$distance)

nrow(HORs.same[HORs.same$chrB == 5,]) / nrow(HORs.same[HORs.same$chrB == 5,])
length(HORs$index[HORs$chrA == 1]) - nrow(HORs.same[HORs.same$chrB == 1,])

mean(HORs.same$distance[HORs.same$chrA == 1])


#gap analysis

gap.100to10000 = data.frame(chromosome1 = vector(mode = "numeric", length = 0), chromosome2 = vector(mode = "numeric", length = 0), start = vector(mode = "numeric", length = 0), end = vector(mode = "numeric", length = 0))

for(i in 2 : nrow(cen180))
{
  print(i)
  if(((cen180$start[i] - cen180$end[i-1]) > 1000 ))
  {
    temp = data.frame(chromosome1 = cen180$chromosome[i], chromosome2 = cen180$chromosome[i-1], start = cen180$start[i], end = cen180$end[i-1])
    gap.100to10000 = rbind(gap.100to10000, temp)
  }
}

gap.100to10000.same = gap.100to10000[gap.100to10000$chromosome1 == gap.100to10000$chromosome2,]




#################


####telo


TELOplus = tolower("TTTAGGG")
TELOminus = tolower("CCCTAAA")

class.names = c("TELOplus", "TELOminus")

class.sequences = c(TELOplus, TELOminus)

allRepeats = data.frame()
allRepeats$contig = vector(mode = "numeric")
allRepeats$class = vector(mode = "numeric")
allRepeats$start = vector(mode = "numeric")
allRepeats$end = vector(mode = "numeric")
allRepeats$sequence = vector(mode = "character")
allRepeats$strand = vector(mode = "character")

a = read.fasta(file = "C:/Users/wlodz/Desktop/Chr1cen4mb.fasta", as.string = TRUE)


for(i in 1 : 1)
{
  time = Sys.time()
  print(paste("fasta file ", i, "/7", sep = ""))
  for(j in 1 : 2)#length(class.names))
  {
    print(paste("repeat class: ", class.names[j], sep = ""))
    maxMis = 0
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

allRepeats$start = allRepeats$start + 13999999
allRepeats$end = allRepeats$end + 13999999

for( i in 1 : length( allRepeats$x ))
{
  print(i)
  if(allRepeats$class[i] == "TELOminus")
  {
    allRepeats$seqSameOrientation[i] = revCompString(allRepeats$seqSameOrientation[i])
  }
}

plot(allRepeats$start, allRepeats$length)

iii = nrow(allRepeats)

while(iii > 1)
{
  print(iii)
  overlap = allRepeats$start[iii] - allRepeats$end[iii - 1] - 1
  
  if(overlap < -3 & overlap > -500) # if overlap is big, the next seq will be removed completely
  {
    allRepeats = allRepeats[-c(iii),] 
  }
  iii = iii - 1
}
allRepeats = allRepeats[,-c(1,5)]
allRepeats$contig = "chr1"

write.csv(x = allRepeats, file = "C:/Users/wlodz/Desktop/TELOchr1repeat.csv", row.names = FALSE)











################# check Ians counts
racon.cen180$seq2 = ""
  
for(i in 1 : nrow(racon.cen180))
{
  print(i)
  racon.cen180$seq2[i] = substr(a[[racon.cen180$chromosome[i]]][1], start = racon.cen180$start[i], racon.cen180$end[i])
}

for(i in 1 : nrow(racon.cen180))
{
  print(i)
  if(cen180$strand.1[i] == "cen180minus")
  {
    cen180$seq2[i] = revCompString(cen180$seq2[i])
  }
}
cen180$seqSameOrientation = cen180$seq2

uneven1 = cen180[cen180$seq2 != cen180b$seqSameOrientation,]

uneven2 = cen180b[cen180$seq2 != cen180b$seqSameOrientation,]

uneven1$len1 = nchar(uneven1$seq2)
uneven1$len2 = nchar(uneven2$seqSameOrientation)


############### extract TSD flanking sequence from athilla annotations and export to align
setwd("C:/Users/wlodz/Desktop/col0610/")

athilas.col0610 = read.csv("Athila.Table.S4_revision.csv")
o = 999
for(i in 1 : nrow(athilas.col0610))
{
  print(i)
  if(athilas.col0610$TSD[i] == "yes")
  {
    write.fasta(sequences = substr(a[[as.numeric(substr(athilas.col0610$ï..Chromosome[i],4,4))]][1], start = athilas.col0610$Start[i] - o, stop = athilas.col0610$Start[i]), 
                names = paste(i, "plus", athilas.col0610$ï..Chromosome[i], athilas.col0610$Class[i], sep = "_"), 
                file.out = "athila.1000.borders.containTSD.fasta", 
                open = "a", 
                as.string = TRUE)
    write.fasta(sequences = revCompString(substr(a[[as.numeric(substr(athilas.col0610$ï..Chromosome[i],4,4))]][1], start = athilas.col0610$End[i], stop = athilas.col0610$End[i] + o)), 
                names = paste(i, "minus",athilas.col0610$ï..Chromosome[i], athilas.col0610$Class[i], sep = "_"), 
                file.out = "athila.1000.borders.containTSD.fasta", 
                open = "a", 
                as.string = TRUE)
  }
}

for(i in 1 : nrow(athilas.col0610))
{
  print(i)
  write.fasta(sequences = substr(a[[as.numeric(substr(athilas.col0610$ï..Chromosome[i],4,4))]][1], start = athilas.col0610$Start[i] - o, stop = athilas.col0610$Start[i]), 
              names = paste(i, "plus", athilas.col0610$ï..Chromosome[i], athilas.col0610$Class[i], sep = "_"), 
              file.out = "athila.1000.borders.all.fasta", 
              open = "a", 
              as.string = TRUE)
  write.fasta(sequences = revCompString(substr(a[[as.numeric(substr(athilas.col0610$ï..Chromosome[i],4,4))]][1], start = athilas.col0610$End[i], stop = athilas.col0610$End[i] + o)), 
              names = paste(i, "minus", athilas.col0610$ï..Chromosome[i], athilas.col0610$Class[i], sep = "_"), 
              file.out = "athila.1000.borders.all.fasta", 
              open = "a", 
              as.string = TRUE)
}






