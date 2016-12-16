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
boxplot(lengths[[1]],lengths[[3]],main="Length (BP)",border=c("#3288bd","#f46d43"),pch=16)
dev.off()
png("enhancer_gc_pc.png")
boxplot(gc.pc[[1]],gc.pc[[3]],Main="GC%",border=c("#3288bd","#f46d43"),pch=16)
abline(h=global.gc.pc,col="darkgrey",lty=2,lwd=2)
dev.off()
png("enhancer_cg_oe.png")
boxplot(cg.oe[[1]],cg.oe[[3]],main="CG OE",border=c("#3288bd","#f46d43"),pch=16)
abline(h=global.oe,col="darkgrey",lty=2,lwd=2)
dev.off()
############################################################################
#
## Manuscript figure 2
#
############################################################################
#
#
#
############################################################################
#
## Manuscript figure 3
#
############################################################################
library(vioplot)
setwd("X:/Current lab members/Emma Bell/Data/Collaborators/Alice Jouneau/RNA-seq")
expr <- read.table("expression_of_SE_associated_genes_vitro_vivo.txt",sep="\t",head=T,comment.char="",quote="",stringsAsFactors=F)
# Each sample is provided in triplicate
# Calculate the geometric mean and standard deviation across each sample
geomean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}
# Samples: 
# 	1. ESC.2i.BU
#	2. EpiLC.BU
#	3. ESC.2i.MA
#	4. ESC.ser.MA
#	5. EpiSC
#	6. ICM
#	7. ERSE
#	8. APE
expr.comb <- array(NA,dim=c(nrow(expr),10))
expr.comb <- data.frame(expr.comb)
colnames(expr.comb) <- c(colnames(expr)[1:2],"ESC.2i.BU","EpiLC.BU","ESC.2i.MA","ESC.ser.MA","EpiSC","ICM","ERSE","APE")
expr.comb[,1:2] <- expr[,1:2]
reps <- list(3:5,6:8,9:12,13:15,16:18,19:21,22:24,25:27)
for(i in 1:8){
y <- i+2
for(j in 1:nrow(expr)){
expr.comb[j,y] <- mean(as.numeric(expr[j,as.numeric(reps[[i]])],na.rm=T))
}
}
# Offset by 1 to allow log transformation
expr.comb[,3:10] <- expr.comb[,3:10]+1
expr.comb[which(expr.comb[,1]=="Mixed"|expr.comb[,1]=="Hypomethylated"),1] <- "Maintained"
boxplot(log10(expr.comb[which(expr.comb[,1]=="Hypermethylated"),5:7]))
meth <-  na.omit(log10(expr.comb[which(expr.comb[,1]=="Hypermethylated"),5:7]))
main <-  na.omit(log10(expr.comb[which(expr.comb[,1]=="Maintained"),5:7]))
dev.new(width=4,height=6)
vioplot(meth[,1],meth[,2],meth[,3],col="white",rectCol="darkgrey",border="magenta")
vioplot(main[,1],main[,2],main[,3],col="white",rectCol="darkgrey",border=rgb(0,1,0))



