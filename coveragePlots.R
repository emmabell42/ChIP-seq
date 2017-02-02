nice R
require(ShortRead)
library(GenomicRanges)
library(GenomicAlignments)
#
## function to plot coverage across set of defined ranges
#
plotCoverage <- function(cov,regions,precision=1000,normalize.widths=TRUE,col="black",lty=1,lwd=1,plotfile=NULL){

	trackVal <- array(0,dim=c(length(regions),precision))
	cat(paste(nrow(trackVal),"sets of regions provided\n"))
	if(length(col)<nrow(trackVal)) col <- rep(col,nrow(trackVal))
	if(length(lty)<nrow(trackVal)) lty <- rep(lty,nrow(trackVal))
	if(length(lwd)<nrow(trackVal)) lwd <- rep(lwd,nrow(trackVal))
	if(!normalize.widths) refMid <- ceiling(precision/2)
	for(i in 1:length(regions)){
		theseSpaces <- intersect(as.character(unique(regions[[i]]@seqnames)),names(cov))
		cat(paste("region set",i,"has",length(theseSpaces),"spaces \n"))
		nrangesByChr <- rep(NA,length(theseSpaces))
		for(k in 1:length(theseSpaces)){
			cat(paste("\t calculating average coverage for",sum(regions[[i]]@seqnames==theseSpaces[k]),"regions on space",theseSpaces[k],"\n"))
			chrCov <- Views(cov[[theseSpaces[k]]],ranges(regions[[i]][regions[[i]]@seqnames==theseSpaces[k]]))
			nrangesByChr[k] <- length(chrCov)
			for(j in 1:length(chrCov)){
				thisCov <- chrCov[[j]]
				if(normalize.widths){
					if(length(thisCov)<precision) thisCov <- Rle(rep(as.numeric(thisCov),each=floor(precision/length(thisCov))))
					indices <- round((0:precision)*(length(thisCov)/precision))
					starts <- indices[1:(length(indices)-1)]+1
					ends <- indices[2:length(indices)]
					newvals <- mean(Views(thisCov,start=starts,end=ends))
					newvals[is.na(newvals)] <- 0
					trackVal[i,] <- trackVal[i,] + newvals
				}
				if(!normalize.widths){
					thisMid <- ceiling(length(thisCov)/2)
					refStart <- max(c(0,refMid-thisMid))+1
					refEnd <- min(c(precision,refMid+length(thisCov)-thisMid))
					thisStart <- max(c(0,thisMid-refMid))+1
					thisEnd <- min(c(length(thisCov),thisMid+precision-refMid))
					cat(paste("filling ref from",refStart,"to",refEnd,"with Cov from",thisStart,"to",thisEnd,"midpoints:",refMid,thisMid,"precision:",precision,"\n"))
					newvals <- as.numeric(thisCov)[thisStart:thisEnd]
					newvals[is.na(newvals)] <- 0
					trackVal[i,refStart:refEnd] <- trackVal[i,refStart:refEnd] + newvals
				}
			}
		}
		trackVal[i,] <- trackVal[i,]/sum(nrangesByChr)
	}

	cat("making plot and returning output \n")
	if(!is.null(plotfile)) png(plotfile)
	cat(paste("range of values:",range(trackVal),"\n"))
	plot(trackVal[1,],lty=lty[1],col=col[1],lwd=lwd[1],xlab="relative position",ylab="averaged coverage",type="l",ylim=range(trackVal,na.rm=TRUE))
	if(length(regions)>1){
		for(i in 2:length(regions)){
			points(trackVal[i,],lty=lty[i],col=col[i],lwd=lwd[i],type="l")
		}
	}
	if(!is.null(plotfile)) dev.off()
	trackVal
}
#
## \n
#
bamFiles <- list.files(recursive=T)[grep(".subset.bam",list.files(recursive=T))]
toExclude <- c(".txt",".png","RNA",".pdf","Otx2ko","MEF","C2C12")
bamFiles[grep(paste(toExclude,collapse="|"),bamFiles)] <- NA
bamFiles <- na.omit(bamFiles)
covObjects <- c()
for(i in 1:length(bamFiles)){
last <- length(strsplit(bamFiles[i],"/")[[1]])
covObjects[i] <- strsplit(bamFiles[i],"/")[[1]][last]
}
covObjects <- gsub(".sorted.subset.bam",".cov",covObjects)
covObjects <- gsub(".subset.bam",".cov",covObjects)
for(i in 1:length(bamFiles)){
cat(i,"Reading and removing duplicates from",bamFiles[i],"at",as.character(Sys.time()),"\n",sep=" ")
tmp <- readGAlignments(bamFiles[i], param=ScanBamParam(what="qname"))
tmp <- tmp[!duplicated(tmp@start)]
cat(i,"Coverting",bamFiles[i],"to GRanges and Coverage object at",as.character(Sys.time()),"\n",sep=" ")
tmp.ranges <- as(tmp, "GRanges")
tmp.cov <- coverage(tmp.ranges)
cat(i,"Normalising",bamFiles[i],"to sequencing depth at",as.character(Sys.time()),"\n",sep=" ")
nReads <- length(tmp.ranges)
normFactor <- 10^6/nReads
normCoverage <- tmp.cov*normFactor
assign(covObjects[i],normCoverage)
}
inputs <- cbind(ip=covObjects,inputs="CTRL")
toExclude <- c("WCE","Control","input","Input","FAIRE","ATAC")
controls <- unique(covObjects[grep(paste(toExclude,collapse="|"),covObjects)]) 
inputs[grep(paste(toExclude,collapse="|"),inputs[,1])] <- NA
inputs <- na.omit(inputs)
inputs[grep("Bruce4",inputs[,1]),2] <- "ES-Bruce4-Control-mouse.cov"
inputs[grep("PennState",inputs[,1]),2] <- "ES-E14-Control-mouse-PennState.cov"
inputs[grep("Stanford",inputs[,1]),2] <- "ES-E14-Control-mouse-Stanford.cov"
inputs[grep("EsrrbNULL",inputs[,1]),2] <- "EsrrbNULL_input.cov"
inputs[grep("EsrrbWT",inputs[,1]),2] <- "EsrrbWT_input.cov"
inputs[grep("GSE49848",inputs[,1]),2] <- "GSE49848_Input.cov"
inputs[grep("GSE27841",inputs[,1]),2] <- "GSE27841_mESC_WCE.cov"
inputs[grep("GSE27841_H3K4me1",inputs[,1]),2] <- "GSE27841_mESC_WCE_for_H3K4me1.cov"
inputs[grep("GSE37262",inputs[,1]),2] <- "GSE37262_mESC_Input.cov"
inputs[grep("GSE39610",inputs[,1]),2] <- "GSE39610_MmES_Input.cov"
inputs[grep("GSE61188",inputs[,1]),2] <- "GSE61188_ChIP_WCE.cov"
inputs[grep("GSE56098",inputs[,1]),2] <- "GSE56098_Input_ESC.cov"
inputs[grep("EpiLC_minus",inputs[,1]),2] <- "GSE56098_Input_EpiLC_minusActivinA.cov"
inputs[grep("EpiLC_plus",inputs[,1]),2] <- "GSE56098_Input_EpiLC_plusActivinA.cov"
inputs[grep("GSE22562",inputs[,1]),2] <- "GSE22562_mESC_WCE.cov"
inputs[grep("GSE24164",inputs[,1]),2] <- "GSE24164_mESC_WCE.cov"
inputs[grep("GSE44286",inputs[,1]),2] <- "GSE44286_mESC_WCE.cov"
stanford <- c("MAFK","CHD2","HCFC1","ZNF384","ZC3H11A","Stanford")
inputs[c(grep(paste(stanford,collapse="|"),inputs[,1])),2] <-  "ES-E14-Control-mouse-Stanford.cov"
inputs[grep("CTRL",inputs[,2]),] <- "ES-E14-Control-mouse-PennState.cov"
for(i in 1:nrow(inputs)){
ip <- get(inputs[i,1])
ip <- ip+1
input <- get(inputs[i,2])
input <- input+1
ip.norm <- ip/input
assign(gsub(".cov",".norm.cov",inputs[i,1]),ip.norm)
}
#
## Read in regions
#
setwd("/data/emmabell42/resources")
bedFiles <- c("DMRandPU_col.bed","DMRregions.txt","ESC_Super_enhancers.bed","protectedregions.txt","SE_meth.bed","TE_coords.bed","ESC_Super_Enhancers_plus20_sorted.bed","ESC_Typical_Enhancers_plus20_sorted.bed")
bedToName <- c("unmeth","hyper","se","hypo","meth","te","sePlus20","tePlus20")
bed.gr <- paste0(bedToName,".gr")
for(i in 1:length(bedFiles)){
bed <- read.table(bedFiles[i],head=F,sep="\t",stringsAsFactors=T)
colnames(bed)[1] <- "chr"
colnames(bed)[2] <- "start"
colnames(bed)[3] <- "end"
assign(bedToName[i],bed)
gr <- as(bed,"GRanges")
assign(bed.gr[i],gr)
}
#
## Calculate coverage
#
toPlot <- c(ls()[grep(".norm.cov",ls())],
ls()[grep("FAIRE",ls())],
ls()[grep("ATAC",ls())])
ChIP.coverage <- array(NA,dim=c(8,1000))
rownames(ChIP.coverage) <- bed.gr
colnames(ChIP.coverage) <- 1:1000
for(i in 1:length(toPlot)){
baseName <- gsub(".norm.cov","",toPlot[i])
rName <- paste0(baseName,".ipnorm.aveCov")
fileName <- paste0(rName,".ipnorm.txt")
ChIP.coverage <- plotCoverage(cov=get(toPlot[i]),regions=list(get(bed.gr[1]),get(bed.gr[2]),get(bed.gr[3]),get(bed.gr[4]),get(bed.gr[5]),get(bed.gr[6]),get(bed.gr[7]),get(bed.gr[8])),normalize.widths=TRUE,precision=1000)
rownames(ChIP.coverage) <- bed.gr
colnames(ChIP.coverage) <- 1:1000
assign(rName,ChIP.coverage)
write.table(ChIP.coverage,fileName,sep="\t",quote=F)
}
#
## Calculate loess line
#
comparisons <- list(SeVsTe=NA,UnmethVsMeth=NA,HypoVsHyperVsMeth=NA)
loessLines <- list(SeVsTe=NA,UnmethVsMeth=NA,HypoVsHyperVsMeth=NA)
x.axis <- 1:1000
for(i in 1:length(toPlot)){
calcLoess <- get(toPlot[i])
comparisons <- list(seVsTe=toPlot[c(3,6),],unmethVsMeth=toPlot[c(1,5),],hypoVsHyperVsMeth=toPlot[c(4,2,5),])
	for(j in 1:3){
	toCompare <- comparisons[[j]]
	loessLines[[j]] <- t(apply(calcLoess,1,function(x) predict(loess(as.numeric(x)~x.axis))))
assign(paste0(toPlot[i],".loess"),loessLines)
}
#
## Plot coverage of:
		#SE vs TE
		#Unmeth vs Meth
		#Hypo vs Hyper vs Meth
#
for(i in 1:length(toPlot)){
toName <- gsub(".cov","",toPlot[i])
covPlot <- get(toPlot[i])
	#SE vs TE
png(paste0(toName,"_SEvsTE.png"),width=240,height=480)
plotCoverage(cov=covPlot,regions=list(se.gr,te.gr),col=c("pink","grey"))
title(toName)
dev.off()
	#Unmeth vs Meth
png(paste0(toName,"_UnmethvsMeth.png"),width=240,height=480)
plotCoverage(cov=covPlot,regions=list(unmeth.gr,meth.gr),col=c("lightblue","grey"))
title(toName)
dev.off()
	#Hypo vs Hyper vs Meth
png(paste0(toName,"_HypovsHypervsMeth.png"),width=240,height=480)
plotCoverage(cov=covPlot,regions=list(hypo.gr,hyper.gr,meth.gr),col=c(rgb(161,215,106,maxColorValue=255),rgb(233,163,201,maxColorValue=255),"grey"))
title(toName)
dev.off()
}
