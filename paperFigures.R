############################################################################
#
## Manuscript figure 1
#
############################################################################
#
## A
#
library(vioplot)
meth <- read.table("X:/Current lab members/Emma Bell/Data/Collaborators/Alice Jouneau/DNA methylation/Promoter methylation/ES_prom_super_enhancer_methylation.txt",sep="\t",head=T,comment.char="",quote="")
dev.new(width=400,height=600)
vioplot(na.omit(meth[,1]),na.omit(meth[,2]),na.omit(meth[,3]),names=c("ES-2i","ES-serum","EpiSC"),col="white",border="black",rectCol="darkgrey")
vioplot(na.omit(meth[,4]),na.omit(meth[,5]),na.omit(meth[,6]),names=c("ES-2i","ES-serum","EpiSC"),col="white",border="black",rectCol="darkgrey")
#
## B
#
setwd("X:/Current lab members/Emma Bell/Data/Collaborators/Alice Jouneau/DNA methylation/Percentage DNA methylation by region ES 2i, serum, EpiSC and ProB cell")
toRead <- list.files()[grep("20161115",list.files())]
for(i in 1:length(toRead)){
tmp <- read.table(toRead[i],sep="\t",head=T,comment.char="",quote="")
assign(toRead[i],tmp)
}
dev.new(width=400,height=600)
vioplot(na.omit(get(toRead[6])[,15]),na.omit(get(toRead[1])[,15]),na.omit(get(toRead[4])[,15]),names=c("ES-2i","ES-serum","EpiSC"),col="white",border="black",rectCol="darkgrey")
vioplot(na.omit(get(toRead[5])[,15]),na.omit(get(toRead[7])[,15]),na.omit(get(toRead[3])[,15]),names=c("ES-2i","ES-serum","EpiSC"),col="white",border="black",rectCol="darkgrey")
#
## C
#
## Created on server as the mouse genome packages no longer work my desktop version of R.
#
library(BSgenome.Mmusculus.UCSC.mm9, quietly = TRUE)
library(TxDb.Mmusculus.UCSC.mm9.knownGene, quietly = TRUE)
library(GenomicRanges)
library(stringr)
setwd("/data/emmabell42/resources/")
toRead <- c("ESC_Super_enhancers.bed","TE_coords.bed","Murine proB cell SE locations.txt","C2C12 cell SE locations.txt","Th cell SE locations.txt","Macrophage cell SE locations.txt")
regions <- c("esc.se","esc.te","prob","c2c12","th","macro")
lengths <- as.list(rep(NA,6))
names(lengths) <- regions
cg.content <- lengths
gc.pc <- lengths
cg.oe <- lengths
# Read in bed files, convert to GRanges objects and get sequences
for(i in 1:length(toRead)){
tmp <- read.table(file=toRead[i],sep="\t",comment.char="",quote="",head=F,stringsAsFactors=F)
colnames(tmp) <- c("chr","start","stop")
assign(regions[i],tmp)
# Calculate lengths of enhancers
lengths[[i]] <- tmp[,3]-tmp[,2]
tmp.ranges <- with(tmp, GRanges(chr, IRanges(start=start, end=stop)))
rangesname <- paste0(regions[i],".ranges")
assign(rangesname,tmp.ranges)
seqsname <- paste0(regions[i],".seqs")
tmp.seqs <- getSeq(Mmusculus, seqnames(tmp.ranges), start(tmp.ranges), end(tmp.ranges), as.character = T)
# Calculate CpG, C and G frequency
counts <- cbind(cg=rep(NA,nrow(tmp)),c=rep(NA,nrow(tmp)),g=rep(NA,nrow(tmp)))
for(j in 1:nrow(counts)){
	counts[j,1] <- str_count(tmp.seqs[j], fixed("CG"))
	counts[j,2] <- str_count(tmp.seqs[j], fixed("C"))
	counts[j,3] <- str_count(tmp.seqs[j], fixed("G"))
	}
cg.content[[i]] <- counts	
# calculate GC%
gc.pc[[i]] <- (counts[,2]+counts[,3])/lengths[[i]]
# Calculate observed/expected CpGs with the Saxonov method
cg.oe[[i]] <- (counts[,1]/lengths[[i]])/((gc.pc[[i]]/2)^2)
}
#Global ratio of observed to expected CpGs (0.1908349)
global.gc.pc <- 0.4171 
global.oe <- 0.0083/((0.4171/2)^2)
# Remove the huge outlier in the TE
cg.oe[[2]][which(cg.oe[[2]]>2)] <- NA
png("enhancer_lengths.png")
boxplot(lengths[[1]],lengths[[3]],main="Length")
dev.off()
png("enhancer_gc_pc.png")
boxplot(gc.pc[[1]],gc.pc[[3]],Main="GC%")
abline(h=global.gc.pc,col="darkgrey",lty=2,lwd=2)
dev.off()
png("enhancer_cg_oe.png")
boxplot(cg.oe[[1]],cg.oe[[3]],main="CG OE")
abline(h=global.oe,col="darkgrey",lty=2,lwd=2)
dev.off()
