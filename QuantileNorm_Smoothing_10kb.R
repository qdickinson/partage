###############################################################################
###############################################################################
#
#		Repli-seq normalization, smoothing and prediction at 10kb resolution
#		November-2021
#
###############################################################################
###############################################################################
#
#	This script normalize RT datasets using Quantile Normalization against a 
#	target dataset (a target dataset would be a high ACF RT dataset with nice 
#	data value distribution), then perform loess smoothing and collect data points
#	at 10Kb windows.
#
#	Requirements:
#	1) RT dataset(s)
#	2) A target dataset
#	3) A file with 1Kb windows of the reference genome
#	4) A file with the chromosome positions to conver hg19 to hg38 (if needed)
#	5) R libraries:
#		- Limma =3.54.2
#		- R.utils =2.12.3
#		- preprocessCore =1.60.2
#		- reshape2 =1.4.4
#		*To install R packages follow these instructions:
#		https://www.msi.umn.edu/support/faq/how-can-i-install-r-packages-my-home-directory

##########################################################################################

library("limma")
library("R.utils")
library("reshape2")

options(scipen=999)

##############  I. VARIABLES TO SET		##################################################

# Write the name you want for your datasets project

ProjectName="RT_20000"

#	If you need to convert from hg19 to Hg38 change this to "T"

convert <- F

##########################################################################################
#
#	DO NOT CHANGE  ANYTHING BELOW THIS POINT!
#
##########################################################################################

fname <- '/home/riveramj/shared/TargetDataset_RT.txt'#target dataset in the format ('Chromosome', 'Position', 'Mean')

rowNum <- countLines(fname)
classes <- NULL
classes[1] = "character"
classes[2:3] = "numeric"
target <- read.delim(fname, header=T, nrows=rowNum, comment.char = "", colClasses=classes)

options("scipen"=20)

##########################################################################################
#
#		Quantile Normalization Function
#

rescaleDatasetsQuantile = function(RTdata, target) {   
  library(preprocessCore)
  if(!missing(target)) {
    ad = target
  } else {
    ad = stack(RTdata[,1:ncol(RTdata)])$values	
  }
  am = as.matrix(RTdata[,1:ncol(RTdata)])
  lq = normalize.quantiles.use.target(am, target=ad)
  RTlq = data.frame(lq)
  names(RTlq) = names(RTdata)
  return(RTlq)
}

##########################################################################################
#
#		Chromosome sizes
#
#hg38 size list
sizelist <- NULL
sizelist <- data.frame(Chromosome = "chr1", Start = 0, End = 248956422, stringsAsFactors = F)
sizelist <- rbind(sizelist, c("chr10", 0, 133797422))
sizelist <- rbind(sizelist, c("chr11", 0, 135086622))
sizelist <- rbind(sizelist, c("chr12", 0, 133275309))
sizelist <- rbind(sizelist, c("chr13", 0, 114364328))
sizelist <- rbind(sizelist, c("chr14", 0, 107043718))
sizelist <- rbind(sizelist, c("chr15", 0, 101991189))
sizelist <- rbind(sizelist, c("chr16", 0, 90338345))
sizelist <- rbind(sizelist, c("chr17", 0, 83257441))
sizelist <- rbind(sizelist, c("chr18", 0, 80373285))
sizelist <- rbind(sizelist, c("chr19", 0, 58617616))
sizelist <- rbind(sizelist, c("chr2", 0, 242193529))
sizelist <- rbind(sizelist, c("chr20", 0, 64444167))
sizelist <- rbind(sizelist, c("chr21", 0, 46709983))
sizelist <- rbind(sizelist, c("chr22", 0, 50818468))
sizelist <- rbind(sizelist, c("chr3", 0, 198295559))
sizelist <- rbind(sizelist, c("chr4", 0, 190214555))
sizelist <- rbind(sizelist, c("chr5", 0, 181538259))
sizelist <- rbind(sizelist, c("chr6", 0, 170805979))
sizelist <- rbind(sizelist, c("chr7", 0, 159345973))
sizelist <- rbind(sizelist, c("chr8", 0, 145138636))
sizelist <- rbind(sizelist, c("chr9", 0, 138394717))
sizelist <- rbind(sizelist, c("chrX", 0, 156040895))
sizelist <- rbind(sizelist, c("chrY", 0, 57227415))
sizelist[,1] <- as.character(sizelist[,1])
sizelist[,2] <- as.numeric(sizelist[,2])
sizelist[,3] <- as.numeric(sizelist[,3])

##########################################################################################
#
#		Chromosome gap list (from ENCODE hg38)
#
#hg38 gap list
gaplist <- NULL
gaplist <- data.frame(Chromosome = "chr1", Start = 122026459, End = 124932724, stringsAsFactors = F)
gaplist <- rbind(gaplist, c("chr10", 39686682, 41593521))
gaplist <- rbind(gaplist, c("chr11", 51078348, 54425074))
gaplist <- rbind(gaplist, c("chr12", 34769407, 37185252))
gaplist <- rbind(gaplist, c("chr13", 16000000, 18051248))
gaplist <- rbind(gaplist, c("chr14", 16000000, 18173523))
gaplist <- rbind(gaplist, c("chr15", 17083673, 19725254))
gaplist <- rbind(gaplist, c("chr16", 36311158, 38265669))
gaplist <- rbind(gaplist, c("chr17", 22813679, 26616164))
gaplist <- rbind(gaplist, c("chr18", 15460899, 20861206))
gaplist <- rbind(gaplist, c("chr19", 24498980, 27190874))
gaplist <- rbind(gaplist, c("chr2", 92188145, 94090557))
gaplist <- rbind(gaplist, c("chr20", 26436232, 28728874))
gaplist <- rbind(gaplist, c("chr21", 10864560, 12915808))
gaplist <- rbind(gaplist, c("chr22", 12954788, 15054318))
gaplist <- rbind(gaplist, c("chr3", 90772458, 93655574))
gaplist <- rbind(gaplist, c("chr4", 49712061, 51743951))
gaplist <- rbind(gaplist, c("chr5", 46485900, 50059807))
gaplist <- rbind(gaplist, c("chr6", 58553888, 59829934))
gaplist <- rbind(gaplist, c("chr7", 58169653, 61528020))
gaplist <- rbind(gaplist, c("chr8", 44033744, 45877265))
gaplist <- rbind(gaplist, c("chr9", 45877265, 45518558))
gaplist <- rbind(gaplist, c("chrX", 58605579, 62412542))
gaplist <- rbind(gaplist, c("chrY", 10316944, 10544039))
gaplist[,1] <- as.character(gaplist[,1])
gaplist[,2] <- as.numeric(gaplist[,2])
gaplist[,3] <- as.numeric(gaplist[,3])



readData <- function(fname, rowNum){
  classes <- NULL
  classes[1] = "character"
  classes[2:4] = "numeric"
  return(read.delim(fname, header=T, nrows=rowNum, comment.char = "", colClasses=classes))
}

readConv <- function(fname, rowNum){
  classes <- NULL
  classes[1:2] = "character"
  classes[3:6] = "numeric"
  return(read.delim(fname, header=T, nrows=rowNum, comment.char = "", colClasses=classes))
}

colunames <- NULL
temp <- NULL
dirlist <- data.matrix(dir(path=getwd(), pattern="RT.bg$", all.files=TRUE, full.names=FALSE, ignore.case=TRUE))
count <- nrow(dirlist)
for (i in 1:count){
  tempname <- strsplit(dirlist[i,1], "RT.bedgraph")
  temp <- paste(tempname[[1]][1], "qNorm","Smoothed",sep="")
  colunames <- cbind(colunames, temp)
}
colu <- NULL
temp <- NULL
for (i in 1: count){
  cat("Reading Data File", i, "\n"); flush.console()
  rowNum <- countLines(dirlist[i,1])
  temp <- list(readData(dirlist[i,1], rowNum))
  colu <- c(colu, temp)
}

ordering <- NULL
ordering <- colu

#=============
#Convert Build
#=============

if (convert){
  thisconv <- "/100718_HG18_WG_CGH_v3.1_HX3_conv_hg19.txt"
  convNum <- countLines(thisconv)
  convert <- readConv(thisconv, convNum)
  
  ordering <- NULL
  for (i in 1:count){
    cat("Converting Data File", i, "\n"); flush.console()
    RT <- data.frame(Chromosome = colu[[i]]$Chromosome, Position = colu[[i]]$Position, Scaled = colu[[i]]$Scaled, stringsAsFactors = F)
    
    placeholder <- RT
    temp <- paste(placeholder$Chromosome, "_", placeholder$Position, sep="")
    placeholder$Score <- temp
    contrast <- convert
    temp <- paste(contrast$Chromosome, "_", contrast$Start, sep="")
    contrast$Score <- temp
    compare1 <- as.matrix(placeholder$Score)
    compare2 <- as.matrix(contrast$Score)
    temp <- compare1[,1] %in% compare2[,1]
    placeholder$keep <- temp
    placeholder <- placeholder[(placeholder$keep == T),]
    placeholder[(ncol(placeholder) - 1):ncol(placeholder)] <- list(NULL)
    placeholder$convStart <- contrast[,5]
    placeholder$convEnd <- contrast[,6]
    RT <- data.frame(Chromosome = placeholder$Chromosome, Position = placeholder$convStart, End = placeholder$convEnd, Scaled = placeholder$Scaled, stringsAsFactors = F)
    temp <- list(RT)
    ordering <- c(ordering, temp)
  }
}

#==========
#File containing predetermined positions

fname <- '/home/riveramj/shared/Predict_1Kb_windows_hg38.txt'
rowNum <- 3088283 #countLines(fname)
classes <- NULL
classes[1] = "character"
classes[2] = "numeric"
predcol <- read.delim(fname, header=T, nrows=rowNum, comment.char = "", colClasses=classes)

#==========

currentDate <- Sys.Date()
temp_name=paste(ProjectName,"_","qNorm","_","Smoothed","_10Kb",currentDate,".txt",sep="")

predorder <- NULL
for (i in 1:count){
  RT <- ordering[[i]]
  colnames(RT)<-c('Chromosome','Position','End','Scaled')
  
  RT$Scaled <- RT$Scaled * 1.59 / IQR(RT$Scaled)
  
  entryi <- data.frame(as.numeric(RT$Scaled))
  entryt <- as.numeric(target$Mean)
  rescaled <- rescaleDatasetsQuantile(entryi, entryt)
  RT <- data.frame(Chromosome = RT$Chromosome, Position = RT$Position, Scaled = as.numeric(rescaled[,1]), stringsAsFactors = F)
  temp_name1=paste(ProjectName,"_","qNorm","_","_10Kb",currentDate,".txt",sep="")
  write.table(RT, file = temp_name1, row.names = F, quote = F, sep = "\t", eol = "\n")
  
  collect <- NULL
  chrs <- unique(RT$Chromosome)
  chrs = chrs[chrs!="chrY"]

  for (chr in chrs){
    cat("Processing Data File", i, "Chromosome", chr, "\n"); flush.console()
    predchr <- predcol[(predcol$Chromosome == chr),]
    thischr <- RT[(RT$Chromosome == chr),]
    chrspan <- 500000 / (max(thischr$Position) - min(thischr$Position))
    thisloess <- loess(thischr$Scaled ~ thischr$Position, span = chrspan)
    
    thispred <- predict(thisloess, predchr$Position)
    thispred[thispred <= -8 | thispred >= 8] <- NA
    predchr$Predict <- thispred
    
    thisgap <- gaplist[(gaplist$Chromosome == chr),]
    predchr$Predict <- ifelse(((predchr$Position >= thisgap$Start) & (predchr$Position <= thisgap$End)), NA, predchr$Predict)
    
    thistest <- NULL
    thiscol <- NULL
    a <- 1
    b <- 10
    while (b <= nrow(predchr)){
      temp <- predchr[a:b,]$Predict
      thistest <- rbind(thistest, temp)
      thiscol <- rbind(thiscol, predchr[a,2])
      a <- (a + 10)
      b <- (b + 10)
    }
    temp <- rowMeans(thistest)
    temp <- data.frame(Chromosome = chr, Position = thiscol, Mean = temp)
    collect <- rbind(collect, temp)
  }
  temp <- list(collect)
  predorder <- c(predorder, temp)
}		


RT <- data.frame(Chromosome = predorder[[1]]$Chromosome, Position = predorder[[1]]$Position, stringsAsFactors = F)
for (i in 1:count){
  RT <- cbind(RT, predorder[[i]]$Mean)
  colnames(RT)[i+2] <- colunames[i]
}

for (g in 1:nrow(gaplist)){
  cat("Gap", gaplist[g,]$Chromosome, "\n"); flush.console()
  RT <- RT[!(RT$Chromosome == gaplist[g,]$Chromosome & ((RT$Position >= gaplist[g,]$Start) & (RT$Position <= gaplist[g,]$End))), ]
}


##########################################################################################


RT <- na.omit(RT)

write.table(RT, file = temp_name, row.names = F, quote = F, sep = "\t", eol = "\n")

RTlong<-melt(RT,id.vars=c("Chromosome", "Position"), variable.name = "Dataset", value.name = "RT")

spt1<-
  lapply(
    split(RTlong,RTlong$Dataset),
    function(x){
      data.frame(Chromosome = x$Chromosome, Start = x$Position, End = x$Position+10000, RT = x$RT)
    }
  )

lapply(names(spt1), function(x){write.table(spt1[[x]], file = paste(substr(x,1,nchar(x)),".bedgraph",sep = ""),row.names=F,col.names = F,quote=F,sep="\t")})

