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
vioplot(na.omit(meth[,1]),na.omit(meth[,2]),na.omit(meth[,3]),names=c("ES-2i","ES-serum","EpiSC"),col="grey",border="black",rectCol="darkgrey")
vioplot(na.omit(meth[,4]),na.omit(meth[,5]),na.omit(meth[,6]),names=c("ES-2i","ES-serum","EpiSC"),col="grey",border="black",rectCol="darkgrey")
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
vioplot(na.omit(get(toRead[6])[,15]),na.omit(get(toRead[1])[,15]),na.omit(get(toRead[4])[,15]),names=c("ES-2i","ES-serum","EpiSC"),col="grey",border="black",rectCol="darkgrey")
vioplot(na.omit(get(toRead[5])[,15]),na.omit(get(toRead[7])[,15]),na.omit(get(toRead[3])[,15]),names=c("ES-2i","ES-serum","EpiSC"),col="grey",border="black",rectCol="darkgrey")
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
# c("#3288bd","#f46d43")
png("enhancer_lengths.png")
boxplot(lengths[[1]],lengths[[2]],lengths[[3]],main="Length (BP)",col="grey",pch=16)
dev.off()
png("enhancer_gc_pc.png")
boxplot(gc.pc[[1]],gc.pc[[2]],gc.pc[[3]],main="GC%",col="grey",pch=16)
abline(h=global.gc.pc,col="darkgrey",lty=2,lwd=2)
dev.off()
png("enhancer_cg_oe.png")
boxplot(cg.oe[[1]],cg.oe[[2]],cg.oe[[3]],main="CG OE",col="grey",pch=16)
abline(h=global.oe,col="darkgrey",lty=2,lwd=2)
dev.off()
############################################################################
#
## Manuscript figure 2
#
############################################################################
#
## C
#
library(GenomicRanges)
subreg <- read.table("SuperEnhancer_subregion_methstatus.txt",sep="\t",head=T,comment.char="",quote="")
subregDMR <- subreg[which(subreg$Status=="DMR"),]
subregProtected <- subreg[which(subreg$Status=="Protected"),]
subregDMR.ranges <- with(subregDMR, GRanges(SE_chr, IRanges(start=Subregion_start, end=Subregion_end)))
subregProtected.ranges <- with(subregProtected, GRanges(SE_chr, IRanges(start=Subregion_start, end=Subregion_end)))
subreg.ranges <- with(subreg, GRanges(Subregion_chr, IRanges(start=Subregion_start, end=Subregion_end)))
subregESCunmeth <- rbind(subregDMR,subregProtected)
subregESCunmeth.ranges <- with(subregESCunmeth, GRanges(Subregion_chr, IRanges(start=Subregion_start, end=Subregion_end)))
subregMeth <- read.table("SE_meth.bed",sep="\t",comment.char="",quote="",head=F,col.names=c("chr","starts","ends","name"),stringsAsFactors=F)
subregMeth.ranges <- with(meth, GRanges(chr, IRanges(start=starts, end=ends)))

DMRlengths <- subregDMR$Subregion_end-subregDMR$Subregion_start
unmethlengths <- subregProtected$Subregion_end-subregProtected$Subregion_start
methlengths <- subregMeth[,3]-subregMeth[,2]

#Does CpG density vary across subregions?

library(BSgenome.Mmusculus.UCSC.mm9, quietly = TRUE)
library(TxDb.Mmusculus.UCSC.mm9.knownGene, quietly = TRUE)

subregDMR.seqs <- getSeq(Mmusculus, seqnames(subregDMR.ranges), start(subregDMR.ranges), end(subregDMR.ranges), as.character = T)

library(stringr)
subregDMR.cg.count <- rep(NA,nrow(subregDMR))
subregDMR.c.count <- rep(NA,nrow(subregDMR))
subregDMR.g.count <- rep(NA,nrow(subregDMR))

for(i in 1:nrow(subregDMR)){
subregDMR.cg.count[i] <- str_count(subregDMR.seqs[i], fixed("CG"))
subregDMR.c.count[i] <- str_count(subregDMR.seqs[i], fixed("C"))
subregDMR.g.count[i] <- str_count(subregDMR.seqs[i], fixed("G"))
}

subregDMR.GC <- (subregDMR.c.count+subregDMR.g.count)/DMRlengths
subregDMR.oe <- (subregDMR.cg.count/(DMRlengths))/((subregDMR.GC/2)^2)

#Plot the CpGs/dinucleotide
boxplot(subregProtected.cg.count/(unmethlengths-1),subregDMR.cg.count/(DMRlengths-1),subregMeth.cg.count/(methlengths-1),names=c("PU","DMR","PM"),col="lightblue",notch=T,ylim=c(0,0.22),ylab="Observed frequency of CpGs")
abline(h=0.0083,col="darkgrey",lty=2)

#Plot lengths of subregions
png("subreg_lengths.png",w=240)
boxplot(unmethlengths,DMRlengths,methlengths,notch=T,names=c("PU","DMR","PM"),col=c("green","magenta","grey"),pch=20)
dev.off()

#Plot observed:expected CpGs
png("subreg_cgoe.png",w=240)
boxplot(subregProtected.oe,subregDMR.oe,subregMeth.oe,names=c("PU","DMR","PM"),col=c("green","magenta","grey"),notch=T,ylab="Observed:expected frequency of CpGs",pch=20)
abline(h=(0.0083/((0.4171/2)^2)),col="darkgrey",lty=2)
dev.off()
#
## E
#
meth <- read.table("X:\\Current lab members\\Emma Bell\\Data\\Collaborators\\Alice Jouneau\\DNA methylation\\SuperEnhancelet_methstatus_PU_DMR_PM_ES_EpiSC.txt",head=T,sep="\t",comment.char="",quote="")
library(vioplot)
dataSets <- c("Marks_2i","Marks_serum","Schubeler_serum","Veillard_EpiSC","Surani_EpiLC","Surani_EpiSC")
methtable <- array(NA,dim=c(6*nrow(meth),6))
methtable <- data.frame(methtable)
colnames(methtable) <- c("SE","Region","Type","Length","DataSet","Meth")
methtable[,1] <- rep(meth[,1],6)
methtable[,2] <- rep(meth[,2],6)
methtable[,3] <- rep(meth[,3],6)
methtable[,4] <- rep(meth[,4],6)
methtable[,5] <- rep(dataSets,each=nrow(meth))
methtable[,6] <- as.vector(as.matrix(meth[,c(5,7,9,11,13,15)]))
for(i in 1:length(dataSets)){
toPlot <- methtable[which(methtable[,5]==dataSets[i]),]
toPlot <- na.omit(toPlot)
fileName <- paste0("SElets_",dataSets[i],".png")
png(fileName,w=240)
boxplot(toPlot[which(toPlot[,3]=="PU"),6],toPlot[which(toPlot[,3]=="DM"),6],toPlot[which(toPlot[,3]=="PM"),6],names=c("PU","DM","PM"),col=c("green","magenta","grey"),pch=20)
title(main=dataSets[i])
dev.off()
}
#
## F
#
meth <- read.table("X:/Current lab members/Emma Bell/Data/Collaborators/Alice Jouneau/DNA methylation/RRBS/synthese_by_CpG.txt",head=T,stringsAsFactors=F,comment.char="",quote="")
samples <- colnames(meth)[seq(7,33,2)]
samples <- gsub("_.*","",samples)
icm.mean <- c()
e6.5.mean <- c()
e7.5.mean <- c()
for(i in 1:nrow(meth)){
icm.mean <- c(icm.mean,mean(as.numeric(meth[i,c(7,9,11,13,15)]),na.rm=T))
e6.5.mean <- c(e6.5.mean,mean(as.numeric(meth[i,c(17,19,21,23)]),na.rm=T))
e7.5.mean <- c(e7.5.mean,mean(as.numeric(meth[i,c(25,27,29,31,33)]),na.rm=T))
}
methtable <- cbind(meth[,1:5],rep(c("ICM","Epiblast_E6.5","Epiblast_E7.5"),each=10662),c(icm.mean,e6.5.mean,e7.5.mean))
methtable <- methtable[which(methtable[,5]==1),]
methtable[which(methtable[,3]=="PU"),3] <- 1
methtable[which(methtable[,3]=="DM"),3] <- 2
methtable[which(methtable[,3]=="PM"),3] <- 3
colnames(methtable)[6:7] <- c("Sample","Meth")
png("ICM.png",w=240)
boxplot(methtable[which(methtable$Sample=="ICM"),7]~methtable[which(methtable$Sample=="ICM"),3],col=c("green","magenta","grey"),main="ICM",pch=20)
dev.off()
png("Epiblast_E6.5.png",w=240)
boxplot(methtable[which(methtable$Sample=="Epiblast_E6.5"),7]~methtable[which(methtable$Sample=="ICM"),3],col=c("green","magenta","grey"),main="Epiblast E6.5",pch=20)
dev.off()
png("Epiblast_E7.5.png",w=240)
boxplot(methtable[which(methtable$Sample=="Epiblast_E7.5"),7]~methtable[which(methtable$Sample=="Epiblast_E7.5"),3],col=c("green","magenta","grey"),main="Epiblast E7.5",pch=20)
dev.off()
############################################################################
#
## Manuscript figure 3
#
############################################################################
#
## B and C
#
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
png("violinplot_rnaseq_invitro_meth.png",w=240)
vioplot(meth[,1],meth[,2],meth[,3],col="pink")
dev.off()
png("violinplot_rnaseq_invitro_main.png",w=240)
vioplot(main[,1],main[,2],main[,3],col="lightgreen")
dev.off()
meth <-  na.omit(log10(expr.comb[which(expr.comb[,1]=="Hypermethylated"),8:10]))
main <-  na.omit(log10(expr.comb[which(expr.comb[,1]=="Maintained"),8:10]))
png("violinplot_rnaseq_invivo_meth.png",w=240)
vioplot(meth[,1],meth[,2],meth[,3],col="pink")
dev.off()
png("violinplot_rnaseq_invivo_main.png",w=240)
vioplot(main[,1],main[,2],main[,3],col="lightgreen")
dev.off()
#
## D 
#
vitro <- read.table("//store.ic.ac.uk/IC/fom/surgeryandcancer/epigenetics-and-development/Current lab members/Emma Bell/Data/Collaborators/Alice Jouneau/RNA-seq/In vitro/transfer_207097_files_9baad94a/ES_serum_vs_EpiSC_result.txt",sep="\t",head=T,stringsAsFactors=F,comment.char="",quote="")
vivo.ape <- read.table("//store.ic.ac.uk/IC/fom/surgeryandcancer/epigenetics-and-development/Current lab members/Emma Bell/Data/Collaborators/Alice Jouneau/RNA-seq/In vitro/transfer_207097_files_9baad94a/mm9 - ICM vs APE - annotated.txt",sep="\t",head=T,stringsAsFactors=F,comment.char="",quote="")
vivo.erse <- read.table("//store.ic.ac.uk/IC/fom/surgeryandcancer/epigenetics-and-development/Current lab members/Emma Bell/Data/Collaborators/Alice Jouneau/RNA-seq/In vitro/transfer_207097_files_9baad94a/mm9 - ICM vs ERSE - annotated.txt",sep="\t",head=T,stringsAsFactors=F,comment.char="",quote="")
pathways <- read.table(file="//store.ic.ac.uk/IC/fom/surgeryandcancer/epigenetics-and-development/Current lab members/Emma Bell/Data/Resources/CPDB_pathways_genes_mouse.txt",sep="\t",header=TRUE)
se <- read.table("//store.ic.ac.uk/IC/fom/surgeryandcancer/epigenetics-and-development/Current lab members/Emma Bell/Data/Analysis/Super Enhancers/Annotation/SE_mRNAs.txt",sep="\t",head=T,stringsAsFactors=F,comment.char="",quote="")
vitro$Adjust.pValue[which(vitro$Adjust.pValue==0)] <- 0.001
selets <- read.table("//store.ic.ac.uk/IC/fom/surgeryandcancer/epigenetics-and-development/Current lab members/Emma Bell/Data/Analysis/Super Enhancers/Super Enhancelets/SuperEnhancer_subregion_methstatus.txt",sep="\t",head=T,stringsAsFactors=F,comment.char="",quote="")
selets$SE.type[which(selets$SE.type=="Hypomethylated")] <- "Maintained"
selets$SE.type[which(selets$SE.type=="Mixed")] <- "Maintained"
selets$SE.type[which(selets$SE.type=="Hypermethylated")] <- "Silenced"
maintained <- unique(selets$Associated.gene[which(selets$SE.type=="Maintained")])
silenced <- unique(selets$Associated.gene[which(selets$SE.type=="Silenced")])
silenced[which(silenced=="")] <- NA
silenced[which(silenced=="Hypermethylatedt1")] <- "Dmrt1"
silenced <- na.omit(silenced)
colnames(vitro)[5] <- "adj.pValue"
colnames(vivo.ape)[2] <- "Symbol"
colnames(vivo.ape)[4] <- "adj.pValue"
colnames(vivo.erse)[5] <- "Symbol"
colnames(vivo.erse)[3] <- "adj.pValue"
pathway.genes <- list()
for(i in 1:nrow(pathways)){
        pathway.genes[[i]] <- strsplit(as.character(pathways$hgnc_symbol_ids[i]),split=",")[[1]]
}
#Create results objects
pathway.enrichments <- array(NA,dim=c(nrow(pathways),3))
rownames(pathway.enrichments) <- pathways$pathway
#Run geneSetTest
experiment <- c("vitro","vivo.erse","vivo.ape")
colnames(pathway.enrichments) <- experiment
library(limma)
for(i in 2:3){
toTest <- get(experiment[i])
for(j in 1:nrow(pathways)){
		# default is for t-statistics. if using f-statistics (as in Fisher's meta-analysis statistics), comment out line 20 and uncomment line 21
                pathway.enrichments[j,i] <- geneSetTest(index=which(toTest$Symbol %in% pathway.genes[[j]]),statistics=(-log10(toTest$adj.pValue)),ranks.only=FALSE,type="f",alternative="mixed")
		#pathway.enrichments[j,i] <- geneSetTest(selected=which(rownames(statistics.table) %in% pathway.genes[[j]]),statistics=statistics.table[,i],type="f",ranks.only=FALSE)
}
pathway.enrichments[,i] <- p.adjust(pathway.enrichments[,i],method="BH")
}
se.geneSetTest <- array(NA,dim=c(3,3))
rownames(se.geneSetTest) <- c("SE","Maintained","Silenced")
colnames(se.geneSetTest) <- experiment
for(i in 1:3){
toTest <- get(experiment[i])
se.geneSetTest[1,i] <- geneSetTest(index=which(toTest$Symbol %in% se[,2]),statistics=(-log10(toTest$adj.pValue)),ranks.only=FALSE,type="f",alternative="mixed")
se.geneSetTest[2,i] <- geneSetTest(index=which(toTest$Symbol %in% maintained),statistics=(-log10(toTest$adj.pValue)),ranks.only=FALSE,type="f",alternative="mixed")
se.geneSetTest[3,i] <- geneSetTest(index=which(toTest$Symbol %in% silenced),statistics=(-log10(toTest$adj.pValue)),ranks.only=FALSE,type="f",alternative="mixed")
se.geneSetTest[,i] <- p.adjust(se.geneSetTest[,i],method="BH")
}
#
## E and F
#
cnrq <- read.table("X:/Current lab members/Emma Bell/Lab work/qPCR data/20150826 AJ conversion data/Final/CNRQ.txt",head=T,stringsAsFactors=F,comment.char="",quote="",sep="\t")
sumtab <- read.table("X:/Current lab members/Emma Bell/Lab work/qPCR data/20150826 AJ conversion data/Final/dataToPlot.txt",head=T,stringsAsFactors=F,comment.char="",quote="",sep="\t")
cnrq[grep("R1",cnrq$Clone),"Clone"] <- "c14"
targets <- unique(cnrq$Target)
clones <- unique(cnrq$Clone)

for(i in 1:length(targets)){
target <- targets[i]
toplot <- cnrq[which(cnrq$Target==target),]
average <- sumtab[which(sumtab$Target==target),]
toplot[which(toplot$Day==10),5] <- 8
average[which(average$Day==10),3] <- 8
ylim.2i <- c(0,max(toplot$RNCQ[which(toplot$Condition=="2i")],na.rm=T)*1.1)
ylim.serum <- c(0,max(toplot$RNCQ[which(toplot$Condition=="Serum")],na.rm=T)*1.1)
xlim <-c(0,8)
png(paste0(target,"_2i_RelExp.png"),w=240)
lp <- plot(RNCQ~Day,data=toplot[which(toplot$Condition=="2i" & toplot$Clone==clones[1] & toplot$Day<8),],xlim=xlim,ylim=ylim.2i,type="o",pch=20,col="grey",main=target,names.arg=NULL,xaxt = 'n',yaxt = 'n',xlab="Day",ylab="Relative expression")
grid(nx=NA,ny=NULL)
points(RNCQ~Day,data=toplot[which(toplot$Condition=="2i" & toplot$Clone==clones[2] & toplot$Day<8),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RNCQ~Day,data=toplot[which(toplot$Condition=="2i" & toplot$Clone==clones[3] & toplot$Day<8),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RNCQ~Day,data=toplot[which(toplot$Condition=="2i" & toplot$Clone==clones[4] & toplot$Day<8),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RNCQ~Day,data=toplot[which(toplot$Condition=="2i" & toplot$Clone==clones[5] & toplot$Day<8),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RNCQ~Day,data=toplot[which(toplot$Condition=="2i" & toplot$Clone==clones[6] & toplot$Day<8),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RelExp~Day,data=average[which(average$Condition=="2i" & average$Day<8),],type="o",pch=16,names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RNCQ~Day,data=toplot[which(toplot$Condition=="2i" & toplot$Day==8),],type="p",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RelExp~Day,data=average[which(average$Condition=="2i" & average$Day==8),],type="p",pch=16,names.arg=NULL,xaxt = 'n',yaxt = 'n')
abline(v=7.5,lty=2,col="darkgrey")
axis(1, at=c(0,1,2,3,7,8), labels=c(0,1,2,3,7,"cEpiSC"),tick=T, las=1, cex.axis=1)
axis(2, at=axTicks(2), cex.axis=0.9, las=2)
dev.off()
png(paste0(target,"_Serum_RelExp.png"),w=240)
lp <- plot(RNCQ~Day,data=toplot[which(toplot$Condition=="Serum" & toplot$Clone==clones[1] & toplot$Day<8),],xlim=xlim,ylim=ylim.serum,type="o",pch=20,col="grey",main=target,names.arg=NULL,xaxt = 'n',yaxt = 'n',xlab="Day",ylab="Relative expression")
grid(nx=NA,ny=NULL)
points(RNCQ~Day,data=toplot[which(toplot$Condition=="Serum" & toplot$Clone==clones[2] & toplot$Day<8),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RNCQ~Day,data=toplot[which(toplot$Condition=="Serum" & toplot$Clone==clones[3] & toplot$Day<8),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RNCQ~Day,data=toplot[which(toplot$Condition=="Serum" & toplot$Clone==clones[4] & toplot$Day<8),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RNCQ~Day,data=toplot[which(toplot$Condition=="Serum" & toplot$Clone==clones[5] & toplot$Day<8),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RNCQ~Day,data=toplot[which(toplot$Condition=="Serum" & toplot$Clone==clones[6] & toplot$Day<8),],type="o",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RelExp~Day,data=average[which(average$Condition=="Serum" & average$Day<8),],type="o",pch=16,names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RNCQ~Day,data=toplot[which(toplot$Condition=="Serum" & toplot$Day==8),],type="p",pch=20,col="grey",names.arg=NULL,xaxt = 'n',yaxt = 'n')
points(RelExp~Day,data=average[which(average$Condition=="Serum" & average$Day==8),],type="p",pch=16,names.arg=NULL,xaxt = 'n',yaxt = 'n')
abline(v=7.5,lty=2,col="darkgrey")
axis(1, at=c(0,1,2,3,7,8), labels=c(0,1,2,3,7,"cEpiSC"),tick=T, las=1, cex.axis=1)
axis(2, at=axTicks(2), cex.axis=0.9, las=2)
dev.off()
}

library(vioplot)
ylim <- extendrange(c(min(log2(na.omit(sumtab$RelExp))),log2(max(sumtab$RelExp[which(sumtab$SE.type=="Maintained")])),f=0.05))
png("Conversion_boxplot_2i_maintained.png",w=240)
boxplot(log2(na.omit(sumtab[which(sumtab$Condition=="2i" & sumtab$SE.type=="Maintained" & sumtab$Day==1),"RelExp"])),
log2(na.omit(sumtab[which(sumtab$Condition=="2i" & sumtab$SE.type=="Maintained" & sumtab$Day==2),"RelExp"])),
log2(na.omit(sumtab[which(sumtab$Condition=="2i" & sumtab$SE.type=="Maintained" & sumtab$Day==3),"RelExp"])),
log2(na.omit(sumtab[which(sumtab$Condition=="2i" & sumtab$SE.type=="Maintained" & sumtab$Day==7),"RelExp"])),
log2(na.omit(sumtab[which(sumtab$Condition=="2i" & sumtab$SE.type=="Maintained" & sumtab$Day==10),"RelExp"])),col="lightgreen",ylim=ylim,pch=20)
dev.off()
png("Conversion_boxplot_2i_silenced.png",w=240)
boxplot(log2(na.omit(sumtab[which(sumtab$Condition=="2i" & sumtab$SE.type=="Silenced" & sumtab$Day==1),"RelExp"])),
log2(na.omit(sumtab[which(sumtab$Condition=="2i" & sumtab$SE.type=="Silenced" & sumtab$Day==2),"RelExp"])),
log2(na.omit(sumtab[which(sumtab$Condition=="2i" & sumtab$SE.type=="Silenced" & sumtab$Day==3),"RelExp"])),
log2(na.omit(sumtab[which(sumtab$Condition=="2i" & sumtab$SE.type=="Silenced" & sumtab$Day==7),"RelExp"])),
log2(na.omit(sumtab[which(sumtab$Condition=="2i" & sumtab$SE.type=="Silenced" & sumtab$Day==10),"RelExp"])),col="pink",ylim=ylim,pch=20)
dev.off()
png("Conversion_boxplot_serum_maintained.png",w=240)
boxplot(log2(na.omit(sumtab[which(sumtab$Condition=="Serum" & sumtab$SE.type=="Maintained" & sumtab$Day==1),"RelExp"])),
log2(na.omit(sumtab[which(sumtab$Condition=="Serum" & sumtab$SE.type=="Maintained" & sumtab$Day==2),"RelExp"])),
log2(na.omit(sumtab[which(sumtab$Condition=="Serum" & sumtab$SE.type=="Maintained" & sumtab$Day==3),"RelExp"])),
log2(na.omit(sumtab[which(sumtab$Condition=="Serum" & sumtab$SE.type=="Maintained" & sumtab$Day==7),"RelExp"])),
log2(na.omit(sumtab[which(sumtab$Condition=="Serum" & sumtab$SE.type=="Maintained" & sumtab$Day==10),"RelExp"])),col="lightgreen",ylim=ylim,pch=20)
dev.off()
png("Conversion_boxplot_Serum_silenced.png",w=240)
boxplot(log2(na.omit(sumtab[which(sumtab$Condition=="Serum" & sumtab$SE.type=="Silenced" & sumtab$Day==1),"RelExp"])),
log2(na.omit(sumtab[which(sumtab$Condition=="Serum" & sumtab$SE.type=="Silenced" & sumtab$Day==2),"RelExp"])),
log2(na.omit(sumtab[which(sumtab$Condition=="Serum" & sumtab$SE.type=="Silenced" & sumtab$Day==3),"RelExp"])),
log2(na.omit(sumtab[which(sumtab$Condition=="Serum" & sumtab$SE.type=="Silenced" & sumtab$Day==7),"RelExp"])),
log2(na.omit(sumtab[which(sumtab$Condition=="Serum" & sumtab$SE.type=="Silenced" & sumtab$Day==10),"RelExp"])),col="pink",ylim=ylim,pch=20)
dev.off()
##########################################
#
## Figure 6
#
##########################################
macs2 callpeak -t merged.bam/EsrrbNULL_Med1.bam -c merged.bam/EsrrbWT_Med1.bam -f BAM -g mm -n 20170621_NULL_Med1 -B
macs2 callpeak -t merged.bam/EsrrbWT_Med1.bam -c merged.bam/EsrrbNULL_Med1.bam -f BAM -g mm -n 20170621_WT_Med1 -B
macs2 callpeak -t merged.bam/EsrrbNULL_H3K27ac.bam -c merged.bam/EsrrbWT_H3K27ac.bam -f BAM -g mm -n 20170621_NULL_H3K27ac -B
macs2 callpeak -t merged.bam/EsrrbWT_H3K27ac.bam -c merged.bam/EsrrbNULL_H3K27ac.bam -f BAM -g mm -n 20170621_WT_H3K27ac -B
 
