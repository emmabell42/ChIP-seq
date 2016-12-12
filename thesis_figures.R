#########################################################################
#
## Chapter 3, figure 2: Transcription factor coverage at Typical Enhancers is similar to coverage at non-functional Super Enhancers
#
#########################################################################
#
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
  assign(rangesname,tmp.ranges)png("enhancer_lengths.png")
boxplot(lengths[[1]],lengths[[2]],lengths[[3]],main="Length (BP)",border="darkgrey",pch=16)
dev.off()
png("enhancer_gc_pc.png")
boxplot(gc.pc[[1]],gc.pc[[2]],gc.pc[[3]],Main="GC%",border="darkgrey",pch=16)
abline(h=global.gc.pc,col="darkgrey",lty=2,lwd=2)
dev.off()
png("enhancer_cg_oe.png")
boxplot(cg.oe[[1]],cg.oe[[2]],cg.oe[[3]],main="CG OE",border="darkgrey",pch=16)
abline(h=global.oe,col="darkgrey",lty=2,lwd=2)
dev.off()

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
boxplot(lengths[[1]],lengths[[2]],lengths[[3]],main="Length (BP)",col="grey",pch=16)
dev.off()
png("enhancer_gc_pc.png")
boxplot(gc.pc[[1]],gc.pc[[2]],gc.pc[[3]],Main="GC%",col="grey",pch=16)
abline(h=global.gc.pc,col="darkgrey",lty=2,lwd=2)
dev.off()
png("enhancer_cg_oe.png")
boxplot(cg.oe[[1]],cg.oe[[2]],cg.oe[[3]],main="CG OE",col="grey",pch=16)
abline(h=global.oe,col="darkgrey",lty=2,lwd=2)
dev.off()
#
#
#########################################################################
#
## Chapter 3, figure 5: Transcription factor binding at mESC Super Enhancers clusters at hypomethylated constituent enhancers.
#
#########################################################################
#
#
setwd("X:/Current lab members/Emma Bell/Data/Collaborators/Peter Rugg-Gunn/20160504 Methylation values/")
toRead <- list.files()[1:8]
assignedName <- gsub(".txt","",toRead)
assignedName <- gsub(" ","_",assignedName)
for(i in 1:length(toRead)){
tmp <- read.table(toRead[i],head=T,sep="\t",comment.char="",quote="",stringsAsFactors=F)
assign(assignedName[i],tmp)
}
#
# Check that the probenames (i.e. the CpGs) in each of the unmethylated and methylated files are identical before proceeding.
#
for(i in 1:4){
unmethylated <- get(assignedName[seq(2,8,2)[i]])
methylated <- get(assignedName[seq(1,7,2)[i]])
cat(identical(unmethylated[,1],methylated[,1]),"\n")
} 
#
# The probenames in all four pairs of files are identical.
#
regions <- c("SuperEnhancer","Hyper","Hypo","Meth")
for(i in 1:4){
unmethylated <- get(assignedName[seq(2,8,2)[i]])
methylated <- get(assignedName[seq(1,7,2)[i]])
totalReads <- unmethylated[,13:60]+methylated[,13:60]
assign(paste0(regions[i],"TotalReads"),totalReads)
proportionMeth <- methylated[,13:60]/totalReads
assign(paste0(regions[i],"proportionMeth"),proportionMeth)
}
#
## Only plot reads with greater than 4 reads
#
toIncl <- as.list(rep(NA,48))
for(i in 1:ncol(totalReads)){
toIncl[[i]] <- if(length(which(totalReads[,i]>3))>0) which(totalReads[,i]>3) else NA
}
SuperEnhancerproportionMeth.5 <- SuperEnhancerproportionMeth
HyperproportionMeth.5 <- HyperproportionMeth
HypoproportionMeth.5 <- HypoproportionMeth
MethproportionMeth.5 <- MethproportionMeth
for(i in 1:ncol(SuperEnhancerTotalReads)){
SuperEnhancerproportionMeth.5[which(SuperEnhancerTotalReads[,i]<5),i] <- NA
HyperproportionMeth.5[which(HyperTotalReads[,i]<5),i] <- NA
HypoproportionMeth.5[which(HypoTotalReads[,i]<5),i] <- NA
MethproportionMeth.5[which(MethTotalReads[,i]<5),i] <- NA
}
#
## Combine replicates
#
combinedMethylated <- array(NA,dim=c(nrow(unmethylated),15))
colnames(combinedMethylated) <-  c("E2.5_Early","E2.5_Late","E3.5_ICM","E3.5_TE","E4.5_FLU_ICM",
"E4.5_FLU_TE","E4.5_IMP_ICM","E4.5_IMP_TE","E5.5_EPI","E5.5_ExE","E5.5_VE","E6.5_EPI","E6.5_ExE","E6.5_VE","Oocytes")
rownames(combinedMethylated) <- unmethylated[,1]
combinedMethylated[,1] <- methylated[,13]+methylated[,14]
combinedMethylated[,2] <- methylated[,15]+methylated[,16]
combinedMethylated[,3] <- methylated[,17]+methylated[,18]+methylated[,19]
combinedMethylated[,4] <- methylated[,20]+methylated[,21]
combinedMethylated[,5] <- methylated[,22]+methylated[,23]+methylated[,24]+methylated[,25]
combinedMethylated[,6] <- methylated[,26]+methylated[,27]+methylated[,28]+methylated[,29]
combinedMethylated[,7] <- methylated[,30]+methylated[,31]+methylated[,32]+methylated[,33]
combinedMethylated[,8] <- methylated[,34]+methylated[,35]+methylated[,36]+methylated[,37]
combinedMethylated[,9] <- methylated[,38]+methylated[,39]+methylated[,40]
combinedMethylated[,10] <- methylated[,41]+methylated[,42]+methylated[,43]
combinedMethylated[,11] <- methylated[,44]+methylated[,45]+methylated[,46]
combinedMethylated[,12] <- methylated[,47]
combinedMethylated[,13] <- methylated[,48]+methylated[,49]
combinedMethylated[,14] <- methylated[,50]+methylated[,51]
combinedMethylated[,15] <- methylated[,52]
combinedUnmethylated <- array(NA,dim=c(nrow(unmethylated),15))
colnames(combinedUnmethylated) <-  c("E2.5_Early","E2.5_Late","E3.5_ICM","E3.5_TE","E4.5_FLU_ICM",
"E4.5_FLU_TE","E4.5_IMP_ICM","E4.5_IMP_TE","E5.5_EPI","E5.5_ExE","E5.5_VE","E6.5_EPI","E6.5_ExE","E6.5_VE","Oocytes")
rownames(combinedUnmethylated) <- unmethylated[,1]
combinedUnmethylated[,1] <- unmethylated[,13]+unmethylated[,14]
combinedUnmethylated[,2] <- unmethylated[,15]+unmethylated[,16]
combinedUnmethylated[,3] <- unmethylated[,17]+unmethylated[,18]+unmethylated[,19]
combinedUnmethylated[,4] <- unmethylated[,20]+unmethylated[,21]
combinedUnmethylated[,5] <- unmethylated[,22]+unmethylated[,23]+unmethylated[,24]+unmethylated[,25]
combinedUnmethylated[,6] <- unmethylated[,26]+unmethylated[,27]+unmethylated[,28]+unmethylated[,29]
combinedUnmethylated[,7] <- unmethylated[,30]+unmethylated[,31]+unmethylated[,32]+unmethylated[,33]
combinedUnmethylated[,8] <- unmethylated[,34]+unmethylated[,35]+unmethylated[,36]+unmethylated[,37]
combinedUnmethylated[,9] <- unmethylated[,38]+unmethylated[,39]+unmethylated[,40]
combinedUnmethylated[,10] <- unmethylated[,41]+unmethylated[,42]+unmethylated[,43]
combinedUnmethylated[,11] <- unmethylated[,44]+unmethylated[,45]+unmethylated[,46]
combinedUnmethylated[,12] <- unmethylated[,47]
combinedUnmethylated[,13] <- unmethylated[,48]+unmethylated[,49]
combinedUnmethylated[,14] <- unmethylated[,50]+unmethylated[,51]
combinedUnmethylated[,15] <- unmethylated[,52]
#
## Calculate methylation proportion of the combined replicates
#
regions <- c("SuperEnhancer","Hyper","Hypo","Meth")
for(i in 1:4){
# Read in unmethylated table
unmethylated <- get(assignedName[seq(2,8,2)[i]])
# Combine replicates
combinedUnmethylated <- array(NA,dim=c(nrow(unmethylated),15))
combinedUnmethylated <- array(NA,dim=c(nrow(unmethylated),15))
combinedUnmethylated[,1] <- unmethylated[,13]+unmethylated[,14]
combinedUnmethylated[,2] <- unmethylated[,15]+unmethylated[,16]
combinedUnmethylated[,3] <- unmethylated[,17]+unmethylated[,18]+unmethylated[,19]
combinedUnmethylated[,4] <- unmethylated[,20]+unmethylated[,21]
combinedUnmethylated[,5] <- unmethylated[,22]+unmethylated[,23]+unmethylated[,24]+unmethylated[,25]
combinedUnmethylated[,6] <- unmethylated[,26]+unmethylated[,27]+unmethylated[,28]+unmethylated[,29]
combinedUnmethylated[,7] <- unmethylated[,30]+unmethylated[,31]+unmethylated[,32]+unmethylated[,33]
combinedUnmethylated[,8] <- unmethylated[,34]+unmethylated[,35]+unmethylated[,36]+unmethylated[,37]
combinedUnmethylated[,9] <- unmethylated[,38]+unmethylated[,39]+unmethylated[,40]
combinedUnmethylated[,10] <- unmethylated[,41]+unmethylated[,42]+unmethylated[,43]
combinedUnmethylated[,11] <- unmethylated[,44]+unmethylated[,45]+unmethylated[,46]
combinedUnmethylated[,12] <- unmethylated[,47]
combinedUnmethylated[,13] <- unmethylated[,48]+unmethylated[,49]
combinedUnmethylated[,14] <- unmethylated[,50]+unmethylated[,51]
combinedUnmethylated[,15] <- unmethylated[,52]
# Read in methylated table
methylated <- get(assignedName[seq(1,7,2)[i]])
# Combine replicates
combinedMethylated <- array(NA,dim=c(nrow(methylated),15))
combinedMethylated[,1] <- methylated[,13]+methylated[,14]
combinedMethylated[,2] <- methylated[,15]+methylated[,16]
combinedMethylated[,3] <- methylated[,17]+methylated[,18]+methylated[,19]
combinedMethylated[,4] <- methylated[,20]+methylated[,21]
combinedMethylated[,5] <- methylated[,22]+methylated[,23]+methylated[,24]+methylated[,25]
combinedMethylated[,6] <- methylated[,26]+methylated[,27]+methylated[,28]+methylated[,29]
combinedMethylated[,7] <- methylated[,30]+methylated[,31]+methylated[,32]+methylated[,33]
combinedMethylated[,8] <- methylated[,34]+methylated[,35]+methylated[,36]+methylated[,37]
combinedMethylated[,9] <- methylated[,38]+methylated[,39]+methylated[,40]
combinedMethylated[,10] <- methylated[,41]+methylated[,42]+methylated[,43]
combinedMethylated[,11] <- methylated[,44]+methylated[,45]+methylated[,46]
combinedMethylated[,12] <- methylated[,47]
combinedMethylated[,13] <- methylated[,48]+methylated[,49]
combinedMethylated[,14] <- methylated[,50]+methylated[,51]
combinedMethylated[,15] <- methylated[,52]
combinedTotalReads <- combinedUnmethylated+combinedMethylated
colnames(combinedTotalReads) <-  c("E2.5_Early","E2.5_Late","E3.5_ICM","E3.5_TE","E4.5_FLU_ICM",
"E4.5_FLU_TE","E4.5_IMP_ICM","E4.5_IMP_TE","E5.5_EPI","E5.5_ExE","E5.5_VE","E6.5_EPI","E6.5_ExE","E6.5_VE","Oocytes")
rownames(combinedTotalReads) <- unmethylated[,1]
assign(paste0(regions[i],"CombinedTotalReads"),combinedTotalReads)
combinedProportionMeth <- combinedMethylated/combinedTotalReads
colnames(combinedProportionMeth) <-  c("E2.5_Early","E2.5_Late","E3.5_ICM","E3.5_TE","E4.5_FLU_ICM",
"E4.5_FLU_TE","E4.5_IMP_ICM","E4.5_IMP_TE","E5.5_EPI","E5.5_ExE","E5.5_VE","E6.5_EPI","E6.5_ExE","E6.5_VE","Oocytes")
rownames(combinedProportionMeth) <- unmethylated[,1]
assign(paste0(regions[i],"combinedProportionMeth"),combinedProportionMeth)
}
#
## Remove CpGs with fewer than 5 reads and plot
#
SuperEnhancercombinedProportionMeth.5 <- SuperEnhancercombinedProportionMeth
HypercombinedProportionMeth.5 <- HypercombinedProportionMeth
HypocombinedProportionMeth.5 <- HypocombinedProportionMeth
MethcombinedProportionMeth.5 <- MethcombinedProportionMeth
for(i in 1:ncol(SuperEnhancercombinedProportionMeth)){
SuperEnhancercombinedProportionMeth.5[which(SuperEnhancerCombinedTotalReads[,i]<5),i] <- NA
HypercombinedProportionMeth.5[which(HyperCombinedTotalReads[,i]<5),i] <- NA
HypocombinedProportionMeth.5[which(HypoCombinedTotalReads[,i]<5),i] <- NA
MethcombinedProportionMeth.5[which(MethCombinedTotalReads[,i]<5),i] <- NA
}
#
## Violin plots
#
library(vioplot)
# Pluripotent cell populations
vioplot(c(na.omit(HypocombinedProportionMeth.5[,3]),na.omit(HypercombinedProportionMeth.5[,3])),
na.omit(MethcombinedProportionMeth.5[,3]),
c(na.omit(HypocombinedProportionMeth.5[,7]),na.omit(HypercombinedProportionMeth.5[,7])),
na.omit(MethcombinedProportionMeth.5[,7]),
names=c("ICM E3.5","ICM E3.5","ICM E4.5","ICM E4.5"),col="grey",rectCol="darkgrey",ylim=c(0,1))
#########################################################################
#
## Chapter 3, figure 4: mESC Super Enhancers acquire CpG methylation as the ICM differentiates to epiblast.
#
#########################################################################
#
#
setwd("X:/Current lab members/Emma Bell/Data/Collaborators/Peter Rugg-Gunn/20160504 Methylation values/")
toRead <- list.files()[1:8]
assignedName <- gsub(".txt","",toRead)
assignedName <- gsub(" ","_",assignedName)
for(i in 1:length(toRead)){
tmp <- read.table(toRead[i],head=T,sep="\t",comment.char="",quote="",stringsAsFactors=F)
assign(assignedName[i],tmp)
}
#
# Check that the probenames (i.e. the CpGs) in each of the unmethylated and methylated files are identical before proceeding.
#
for(i in 1:4){
unmethylated <- get(assignedName[seq(2,8,2)[i]])
methylated <- get(assignedName[seq(1,7,2)[i]])
cat(identical(unmethylated[,1],methylated[,1]),"\n")
} 
#
# The probenames in all four pairs of files are identical.
#
regions <- c("SuperEnhancer","Hyper","Hypo","Meth")
for(i in 1:4){
unmethylated <- get(assignedName[seq(2,8,2)[i]])
methylated <- get(assignedName[seq(1,7,2)[i]])
totalReads <- unmethylated[,13:60]+methylated[,13:60]
assign(paste0(regions[i],"TotalReads"),totalReads)
proportionMeth <- methylated[,13:60]/totalReads
assign(paste0(regions[i],"proportionMeth"),proportionMeth)
}
#
## Only plot reads with greater than 4 reads
#
toIncl <- as.list(rep(NA,48))
for(i in 1:ncol(totalReads)){
toIncl[[i]] <- if(length(which(totalReads[,i]>3))>0) which(totalReads[,i]>3) else NA
}
SuperEnhancerproportionMeth.5 <- SuperEnhancerproportionMeth
HyperproportionMeth.5 <- HyperproportionMeth
HypoproportionMeth.5 <- HypoproportionMeth
MethproportionMeth.5 <- MethproportionMeth
for(i in 1:ncol(SuperEnhancerTotalReads)){
SuperEnhancerproportionMeth.5[which(SuperEnhancerTotalReads[,i]<5),i] <- NA
HyperproportionMeth.5[which(HyperTotalReads[,i]<5),i] <- NA
HypoproportionMeth.5[which(HypoTotalReads[,i]<5),i] <- NA
MethproportionMeth.5[which(MethTotalReads[,i]<5),i] <- NA
}
#
## Combine replicates
#
combinedMethylated <- array(NA,dim=c(nrow(unmethylated),15))
colnames(combinedMethylated) <-  c("E2.5_Early","E2.5_Late","E3.5_ICM","E3.5_TE","E4.5_FLU_ICM",
"E4.5_FLU_TE","E4.5_IMP_ICM","E4.5_IMP_TE","E5.5_EPI","E5.5_ExE","E5.5_VE","E6.5_EPI","E6.5_ExE","E6.5_VE","Oocytes")
rownames(combinedMethylated) <- unmethylated[,1]
combinedMethylated[,1] <- methylated[,13]+methylated[,14]
combinedMethylated[,2] <- methylated[,15]+methylated[,16]
combinedMethylated[,3] <- methylated[,17]+methylated[,18]+methylated[,19]
combinedMethylated[,4] <- methylated[,20]+methylated[,21]
combinedMethylated[,5] <- methylated[,22]+methylated[,23]+methylated[,24]+methylated[,25]
combinedMethylated[,6] <- methylated[,26]+methylated[,27]+methylated[,28]+methylated[,29]
combinedMethylated[,7] <- methylated[,30]+methylated[,31]+methylated[,32]+methylated[,33]
combinedMethylated[,8] <- methylated[,34]+methylated[,35]+methylated[,36]+methylated[,37]
combinedMethylated[,9] <- methylated[,38]+methylated[,39]+methylated[,40]
combinedMethylated[,10] <- methylated[,41]+methylated[,42]+methylated[,43]
combinedMethylated[,11] <- methylated[,44]+methylated[,45]+methylated[,46]
combinedMethylated[,12] <- methylated[,47]
combinedMethylated[,13] <- methylated[,48]+methylated[,49]
combinedMethylated[,14] <- methylated[,50]+methylated[,51]
combinedMethylated[,15] <- methylated[,52]
combinedUnmethylated <- array(NA,dim=c(nrow(unmethylated),15))
colnames(combinedUnmethylated) <-  c("E2.5_Early","E2.5_Late","E3.5_ICM","E3.5_TE","E4.5_FLU_ICM",
"E4.5_FLU_TE","E4.5_IMP_ICM","E4.5_IMP_TE","E5.5_EPI","E5.5_ExE","E5.5_VE","E6.5_EPI","E6.5_ExE","E6.5_VE","Oocytes")
rownames(combinedUnmethylated) <- unmethylated[,1]
combinedUnmethylated[,1] <- unmethylated[,13]+unmethylated[,14]
combinedUnmethylated[,2] <- unmethylated[,15]+unmethylated[,16]
combinedUnmethylated[,3] <- unmethylated[,17]+unmethylated[,18]+unmethylated[,19]
combinedUnmethylated[,4] <- unmethylated[,20]+unmethylated[,21]
combinedUnmethylated[,5] <- unmethylated[,22]+unmethylated[,23]+unmethylated[,24]+unmethylated[,25]
combinedUnmethylated[,6] <- unmethylated[,26]+unmethylated[,27]+unmethylated[,28]+unmethylated[,29]
combinedUnmethylated[,7] <- unmethylated[,30]+unmethylated[,31]+unmethylated[,32]+unmethylated[,33]
combinedUnmethylated[,8] <- unmethylated[,34]+unmethylated[,35]+unmethylated[,36]+unmethylated[,37]
combinedUnmethylated[,9] <- unmethylated[,38]+unmethylated[,39]+unmethylated[,40]
combinedUnmethylated[,10] <- unmethylated[,41]+unmethylated[,42]+unmethylated[,43]
combinedUnmethylated[,11] <- unmethylated[,44]+unmethylated[,45]+unmethylated[,46]
combinedUnmethylated[,12] <- unmethylated[,47]
combinedUnmethylated[,13] <- unmethylated[,48]+unmethylated[,49]
combinedUnmethylated[,14] <- unmethylated[,50]+unmethylated[,51]
combinedUnmethylated[,15] <- unmethylated[,52]
#
## Calculate methylation proportion of the combined replicates
#
regions <- c("SuperEnhancer","Hyper","Hypo","Meth")
for(i in 1:4){
# Read in unmethylated table
unmethylated <- get(assignedName[seq(2,8,2)[i]])
# Combine replicates
combinedUnmethylated <- array(NA,dim=c(nrow(unmethylated),15))
combinedUnmethylated <- array(NA,dim=c(nrow(unmethylated),15))
combinedUnmethylated[,1] <- unmethylated[,13]+unmethylated[,14]
combinedUnmethylated[,2] <- unmethylated[,15]+unmethylated[,16]
combinedUnmethylated[,3] <- unmethylated[,17]+unmethylated[,18]+unmethylated[,19]
combinedUnmethylated[,4] <- unmethylated[,20]+unmethylated[,21]
combinedUnmethylated[,5] <- unmethylated[,22]+unmethylated[,23]+unmethylated[,24]+unmethylated[,25]
combinedUnmethylated[,6] <- unmethylated[,26]+unmethylated[,27]+unmethylated[,28]+unmethylated[,29]
combinedUnmethylated[,7] <- unmethylated[,30]+unmethylated[,31]+unmethylated[,32]+unmethylated[,33]
combinedUnmethylated[,8] <- unmethylated[,34]+unmethylated[,35]+unmethylated[,36]+unmethylated[,37]
combinedUnmethylated[,9] <- unmethylated[,38]+unmethylated[,39]+unmethylated[,40]
combinedUnmethylated[,10] <- unmethylated[,41]+unmethylated[,42]+unmethylated[,43]
combinedUnmethylated[,11] <- unmethylated[,44]+unmethylated[,45]+unmethylated[,46]
combinedUnmethylated[,12] <- unmethylated[,47]
combinedUnmethylated[,13] <- unmethylated[,48]+unmethylated[,49]
combinedUnmethylated[,14] <- unmethylated[,50]+unmethylated[,51]
combinedUnmethylated[,15] <- unmethylated[,52]
# Read in methylated table
methylated <- get(assignedName[seq(1,7,2)[i]])
# Combine replicates
combinedMethylated <- array(NA,dim=c(nrow(methylated),15))
combinedMethylated[,1] <- methylated[,13]+methylated[,14]
combinedMethylated[,2] <- methylated[,15]+methylated[,16]
combinedMethylated[,3] <- methylated[,17]+methylated[,18]+methylated[,19]
combinedMethylated[,4] <- methylated[,20]+methylated[,21]
combinedMethylated[,5] <- methylated[,22]+methylated[,23]+methylated[,24]+methylated[,25]
combinedMethylated[,6] <- methylated[,26]+methylated[,27]+methylated[,28]+methylated[,29]
combinedMethylated[,7] <- methylated[,30]+methylated[,31]+methylated[,32]+methylated[,33]
combinedMethylated[,8] <- methylated[,34]+methylated[,35]+methylated[,36]+methylated[,37]
combinedMethylated[,9] <- methylated[,38]+methylated[,39]+methylated[,40]
combinedMethylated[,10] <- methylated[,41]+methylated[,42]+methylated[,43]
combinedMethylated[,11] <- methylated[,44]+methylated[,45]+methylated[,46]
combinedMethylated[,12] <- methylated[,47]
combinedMethylated[,13] <- methylated[,48]+methylated[,49]
combinedMethylated[,14] <- methylated[,50]+methylated[,51]
combinedMethylated[,15] <- methylated[,52]
combinedTotalReads <- combinedUnmethylated+combinedMethylated
colnames(combinedTotalReads) <-  c("E2.5_Early","E2.5_Late","E3.5_ICM","E3.5_TE","E4.5_FLU_ICM",
"E4.5_FLU_TE","E4.5_IMP_ICM","E4.5_IMP_TE","E5.5_EPI","E5.5_ExE","E5.5_VE","E6.5_EPI","E6.5_ExE","E6.5_VE","Oocytes")
rownames(combinedTotalReads) <- unmethylated[,1]
assign(paste0(regions[i],"CombinedTotalReads"),combinedTotalReads)
combinedProportionMeth <- combinedMethylated/combinedTotalReads
colnames(combinedProportionMeth) <-  c("E2.5_Early","E2.5_Late","E3.5_ICM","E3.5_TE","E4.5_FLU_ICM",
"E4.5_FLU_TE","E4.5_IMP_ICM","E4.5_IMP_TE","E5.5_EPI","E5.5_ExE","E5.5_VE","E6.5_EPI","E6.5_ExE","E6.5_VE","Oocytes")
rownames(combinedProportionMeth) <- unmethylated[,1]
assign(paste0(regions[i],"combinedProportionMeth"),combinedProportionMeth)
}
#
## Remove CpGs with fewer than 5 reads and plot
#
SuperEnhancercombinedProportionMeth.5 <- SuperEnhancercombinedProportionMeth
HypercombinedProportionMeth.5 <- HypercombinedProportionMeth
HypocombinedProportionMeth.5 <- HypocombinedProportionMeth
MethcombinedProportionMeth.5 <- MethcombinedProportionMeth
for(i in 1:ncol(SuperEnhancercombinedProportionMeth)){
SuperEnhancercombinedProportionMeth.5[which(SuperEnhancerCombinedTotalReads[,i]<5),i] <- NA
HypercombinedProportionMeth.5[which(HyperCombinedTotalReads[,i]<5),i] <- NA
HypocombinedProportionMeth.5[which(HypoCombinedTotalReads[,i]<5),i] <- NA
MethcombinedProportionMeth.5[which(MethCombinedTotalReads[,i]<5),i] <- NA
}
#
## Violin plots
#
library(vioplot)
# Pluripotent cell populations
par(mar=c(6,3,2,2))
vioplot(na.omit(SuperEnhancercombinedProportionMeth.5[,3]),
na.omit(SuperEnhancercombinedProportionMeth.5[,7]),
na.omit(SuperEnhancercombinedProportionMeth.5[,9]),
na.omit(SuperEnhancercombinedProportionMeth.5[,12]),names=c("ICM E3.5","ICM E4.5","Epiblast E5.5","Epiblast E6.5"),col="grey",ylim=c(0,1),rectCol="darkgrey")
#########################################################################
#
## Chapter 3, figure 16: mESC Super Enhancers acquire CpG methylation as the ICM differentiates to epiblast.
#
#########################################################################
#
#
setwd("X:/Current lab members/Emma Bell/Data/Collaborators/Peter Rugg-Gunn/20160504 Methylation values/")
toRead <- list.files()[1:8]
assignedName <- gsub(".txt","",toRead)
assignedName <- gsub(" ","_",assignedName)
for(i in 1:length(toRead)){
tmp <- read.table(toRead[i],head=T,sep="\t",comment.char="",quote="",stringsAsFactors=F)
assign(assignedName[i],tmp)
}
#
# Check that the probenames (i.e. the CpGs) in each of the unmethylated and methylated files are identical before proceeding.
#
for(i in 1:4){
unmethylated <- get(assignedName[seq(2,8,2)[i]])
methylated <- get(assignedName[seq(1,7,2)[i]])
cat(identical(unmethylated[,1],methylated[,1]),"\n")
} 
#
# The probenames in all four pairs of files are identical.
#
regions <- c("SuperEnhancer","Hyper","Hypo","Meth")
for(i in 1:4){
unmethylated <- get(assignedName[seq(2,8,2)[i]])
methylated <- get(assignedName[seq(1,7,2)[i]])
totalReads <- unmethylated[,13:60]+methylated[,13:60]
assign(paste0(regions[i],"TotalReads"),totalReads)
proportionMeth <- methylated[,13:60]/totalReads
assign(paste0(regions[i],"proportionMeth"),proportionMeth)
}
#
## Only plot reads with greater than 4 reads
#
toIncl <- as.list(rep(NA,48))
for(i in 1:ncol(totalReads)){
toIncl[[i]] <- if(length(which(totalReads[,i]>3))>0) which(totalReads[,i]>3) else NA
}
SuperEnhancerproportionMeth.5 <- SuperEnhancerproportionMeth
HyperproportionMeth.5 <- HyperproportionMeth
HypoproportionMeth.5 <- HypoproportionMeth
MethproportionMeth.5 <- MethproportionMeth
for(i in 1:ncol(SuperEnhancerTotalReads)){
SuperEnhancerproportionMeth.5[which(SuperEnhancerTotalReads[,i]<5),i] <- NA
HyperproportionMeth.5[which(HyperTotalReads[,i]<5),i] <- NA
HypoproportionMeth.5[which(HypoTotalReads[,i]<5),i] <- NA
MethproportionMeth.5[which(MethTotalReads[,i]<5),i] <- NA
}
#
## Combine replicates
#
combinedMethylated <- array(NA,dim=c(nrow(unmethylated),15))
colnames(combinedMethylated) <-  c("E2.5_Early","E2.5_Late","E3.5_ICM","E3.5_TE","E4.5_FLU_ICM",
"E4.5_FLU_TE","E4.5_IMP_ICM","E4.5_IMP_TE","E5.5_EPI","E5.5_ExE","E5.5_VE","E6.5_EPI","E6.5_ExE","E6.5_VE","Oocytes")
rownames(combinedMethylated) <- unmethylated[,1]
combinedMethylated[,1] <- methylated[,13]+methylated[,14]
combinedMethylated[,2] <- methylated[,15]+methylated[,16]
combinedMethylated[,3] <- methylated[,17]+methylated[,18]+methylated[,19]
combinedMethylated[,4] <- methylated[,20]+methylated[,21]
combinedMethylated[,5] <- methylated[,22]+methylated[,23]+methylated[,24]+methylated[,25]
combinedMethylated[,6] <- methylated[,26]+methylated[,27]+methylated[,28]+methylated[,29]
combinedMethylated[,7] <- methylated[,30]+methylated[,31]+methylated[,32]+methylated[,33]
combinedMethylated[,8] <- methylated[,34]+methylated[,35]+methylated[,36]+methylated[,37]
combinedMethylated[,9] <- methylated[,38]+methylated[,39]+methylated[,40]
combinedMethylated[,10] <- methylated[,41]+methylated[,42]+methylated[,43]
combinedMethylated[,11] <- methylated[,44]+methylated[,45]+methylated[,46]
combinedMethylated[,12] <- methylated[,47]
combinedMethylated[,13] <- methylated[,48]+methylated[,49]
combinedMethylated[,14] <- methylated[,50]+methylated[,51]
combinedMethylated[,15] <- methylated[,52]
combinedUnmethylated <- array(NA,dim=c(nrow(unmethylated),15))
colnames(combinedUnmethylated) <-  c("E2.5_Early","E2.5_Late","E3.5_ICM","E3.5_TE","E4.5_FLU_ICM",
"E4.5_FLU_TE","E4.5_IMP_ICM","E4.5_IMP_TE","E5.5_EPI","E5.5_ExE","E5.5_VE","E6.5_EPI","E6.5_ExE","E6.5_VE","Oocytes")
rownames(combinedUnmethylated) <- unmethylated[,1]
combinedUnmethylated[,1] <- unmethylated[,13]+unmethylated[,14]
combinedUnmethylated[,2] <- unmethylated[,15]+unmethylated[,16]
combinedUnmethylated[,3] <- unmethylated[,17]+unmethylated[,18]+unmethylated[,19]
combinedUnmethylated[,4] <- unmethylated[,20]+unmethylated[,21]
combinedUnmethylated[,5] <- unmethylated[,22]+unmethylated[,23]+unmethylated[,24]+unmethylated[,25]
combinedUnmethylated[,6] <- unmethylated[,26]+unmethylated[,27]+unmethylated[,28]+unmethylated[,29]
combinedUnmethylated[,7] <- unmethylated[,30]+unmethylated[,31]+unmethylated[,32]+unmethylated[,33]
combinedUnmethylated[,8] <- unmethylated[,34]+unmethylated[,35]+unmethylated[,36]+unmethylated[,37]
combinedUnmethylated[,9] <- unmethylated[,38]+unmethylated[,39]+unmethylated[,40]
combinedUnmethylated[,10] <- unmethylated[,41]+unmethylated[,42]+unmethylated[,43]
combinedUnmethylated[,11] <- unmethylated[,44]+unmethylated[,45]+unmethylated[,46]
combinedUnmethylated[,12] <- unmethylated[,47]
combinedUnmethylated[,13] <- unmethylated[,48]+unmethylated[,49]
combinedUnmethylated[,14] <- unmethylated[,50]+unmethylated[,51]
combinedUnmethylated[,15] <- unmethylated[,52]
#
## Calculate methylation proportion of the combined replicates
#
regions <- c("SuperEnhancer","Hyper","Hypo","Meth")
for(i in 1:4){
# Read in unmethylated table
unmethylated <- get(assignedName[seq(2,8,2)[i]])
# Combine replicates
combinedUnmethylated <- array(NA,dim=c(nrow(unmethylated),15))
combinedUnmethylated <- array(NA,dim=c(nrow(unmethylated),15))
combinedUnmethylated[,1] <- unmethylated[,13]+unmethylated[,14]
combinedUnmethylated[,2] <- unmethylated[,15]+unmethylated[,16]
combinedUnmethylated[,3] <- unmethylated[,17]+unmethylated[,18]+unmethylated[,19]
combinedUnmethylated[,4] <- unmethylated[,20]+unmethylated[,21]
combinedUnmethylated[,5] <- unmethylated[,22]+unmethylated[,23]+unmethylated[,24]+unmethylated[,25]
combinedUnmethylated[,6] <- unmethylated[,26]+unmethylated[,27]+unmethylated[,28]+unmethylated[,29]
combinedUnmethylated[,7] <- unmethylated[,30]+unmethylated[,31]+unmethylated[,32]+unmethylated[,33]
combinedUnmethylated[,8] <- unmethylated[,34]+unmethylated[,35]+unmethylated[,36]+unmethylated[,37]
combinedUnmethylated[,9] <- unmethylated[,38]+unmethylated[,39]+unmethylated[,40]
combinedUnmethylated[,10] <- unmethylated[,41]+unmethylated[,42]+unmethylated[,43]
combinedUnmethylated[,11] <- unmethylated[,44]+unmethylated[,45]+unmethylated[,46]
combinedUnmethylated[,12] <- unmethylated[,47]
combinedUnmethylated[,13] <- unmethylated[,48]+unmethylated[,49]
combinedUnmethylated[,14] <- unmethylated[,50]+unmethylated[,51]
combinedUnmethylated[,15] <- unmethylated[,52]
# Read in methylated table
methylated <- get(assignedName[seq(1,7,2)[i]])
# Combine replicates
combinedMethylated <- array(NA,dim=c(nrow(methylated),15))
combinedMethylated[,1] <- methylated[,13]+methylated[,14]
combinedMethylated[,2] <- methylated[,15]+methylated[,16]
combinedMethylated[,3] <- methylated[,17]+methylated[,18]+methylated[,19]
combinedMethylated[,4] <- methylated[,20]+methylated[,21]
combinedMethylated[,5] <- methylated[,22]+methylated[,23]+methylated[,24]+methylated[,25]
combinedMethylated[,6] <- methylated[,26]+methylated[,27]+methylated[,28]+methylated[,29]
combinedMethylated[,7] <- methylated[,30]+methylated[,31]+methylated[,32]+methylated[,33]
combinedMethylated[,8] <- methylated[,34]+methylated[,35]+methylated[,36]+methylated[,37]
combinedMethylated[,9] <- methylated[,38]+methylated[,39]+methylated[,40]
combinedMethylated[,10] <- methylated[,41]+methylated[,42]+methylated[,43]
combinedMethylated[,11] <- methylated[,44]+methylated[,45]+methylated[,46]
combinedMethylated[,12] <- methylated[,47]
combinedMethylated[,13] <- methylated[,48]+methylated[,49]
combinedMethylated[,14] <- methylated[,50]+methylated[,51]
combinedMethylated[,15] <- methylated[,52]
combinedTotalReads <- combinedUnmethylated+combinedMethylated
colnames(combinedTotalReads) <-  c("E2.5_Early","E2.5_Late","E3.5_ICM","E3.5_TE","E4.5_FLU_ICM",
"E4.5_FLU_TE","E4.5_IMP_ICM","E4.5_IMP_TE","E5.5_EPI","E5.5_ExE","E5.5_VE","E6.5_EPI","E6.5_ExE","E6.5_VE","Oocytes")
rownames(combinedTotalReads) <- unmethylated[,1]
assign(paste0(regions[i],"CombinedTotalReads"),combinedTotalReads)
combinedProportionMeth <- combinedMethylated/combinedTotalReads
colnames(combinedProportionMeth) <-  c("E2.5_Early","E2.5_Late","E3.5_ICM","E3.5_TE","E4.5_FLU_ICM",
"E4.5_FLU_TE","E4.5_IMP_ICM","E4.5_IMP_TE","E5.5_EPI","E5.5_ExE","E5.5_VE","E6.5_EPI","E6.5_ExE","E6.5_VE","Oocytes")
rownames(combinedProportionMeth) <- unmethylated[,1]
assign(paste0(regions[i],"combinedProportionMeth"),combinedProportionMeth)
}
#
## Remove CpGs with fewer than 5 reads and plot
#
SuperEnhancercombinedProportionMeth.5 <- SuperEnhancercombinedProportionMeth
HypercombinedProportionMeth.5 <- HypercombinedProportionMeth
HypocombinedProportionMeth.5 <- HypocombinedProportionMeth
MethcombinedProportionMeth.5 <- MethcombinedProportionMeth
for(i in 1:ncol(SuperEnhancercombinedProportionMeth)){
SuperEnhancercombinedProportionMeth.5[which(SuperEnhancerCombinedTotalReads[,i]<5),i] <- NA
HypercombinedProportionMeth.5[which(HyperCombinedTotalReads[,i]<5),i] <- NA
HypocombinedProportionMeth.5[which(HypoCombinedTotalReads[,i]<5),i] <- NA
MethcombinedProportionMeth.5[which(MethCombinedTotalReads[,i]<5),i] <- NA
}
#
## Violin plots
#
library(vioplot)
vioplot(na.omit(HypercombinedProportionMeth.5[,3]),
na.omit(HypercombinedProportionMeth.5[,7]),
na.omit(HypercombinedProportionMeth.5[,9]),
na.omit(HypercombinedProportionMeth.5[,12]),names=c("ICM E3.5","ICM E4.5","Epiblast E5.5","Epiblast E6.5"),col=rgb(1,0,1
),ylim=c(0,1),rectCol="darkgrey")
vioplot(na.omit(HypocombinedProportionMeth.5[,3]),
na.omit(HypocombinedProportionMeth.5[,7]),
na.omit(HypocombinedProportionMeth.5[,9]),
na.omit(HypocombinedProportionMeth.5[,12]),names=c("ICM E3.5","ICM E4.5","Epiblast E5.5","Epiblast E6.5"),col=rgb(0,1,0
),ylim=c(0,1),rectCol="darkgrey")
