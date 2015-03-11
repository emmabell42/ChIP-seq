## change this to what you expect/
ideal.readlength <- 50

ofile <- as.character(read.table("file_to_filter.txt")[[1]])
setwd("Fastq")
library(ShortRead)
reads <- readFastq(ofile)
ideal.readlength <- as.numeric(names(tail(sort(table(width(reads[1:100]))),n=1)))
quals1 <- as(FastqQuality(quality(quality(reads[1:min(c(length(reads),20000000))]))),"matrix")
if(length(reads)>20000000){
        quals2 <- as(FastqQuality(quality(quality(reads[c(20000001:length(reads))]))),"matrix")
}
trim <- apply(quals1,MARGIN=1,function(x)min(which(cumsum(as.numeric(x<20))>5)))
if(length(reads)>20000000){
        trim2 <- apply(quals2,MARGIN=1,function(x)min(which(cumsum(as.numeric(x<20))>5)))
        trim <- c(trim,trim2)
}
bad.q <- which(!is.infinite(trim) & trim<(0.8*ideal.readlength))
good.reads <- setdiff(c(1:length(reads)),bad.q)
reads <- reads[good.reads]
trim <- trim[good.reads]
trim[is.infinite(trim)] <- ideal.readlength
reads <- ShortReadQ(sread=subseq(sread(reads),start=rep(1,length(trim)),end=trim),quality=new(Class=class(quality(reads)),quality=subseq(quality(quality(reads)),start=rep(1,length(trim)),end=trim)),id=id(reads))
no_n <- which(alphabetFrequency(sread(reads[1:100]))[,'N']==0)
newfile <- strsplit(ofile,split=".fastq")[[1]][1]
newfile <- paste(newfile,"_filtered.fq",sep="")
writeFastq(reads[no_n],newfile)
