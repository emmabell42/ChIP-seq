#ChIP-seq read coverage heatmaps
all.table <- as.matrix(read.table("ChIPseq_coverage_allfactors_allregions.txt",sep="\t",comment.char="",quote="",head=T))
colnames(all.table) <- 1:1000

esc <- all.table[grep("mESC",rownames(all.table)),]
epi <- all.table[grep("Epi",rownames(all.table)),]

png("esc_dendro.png")
plot(hclust(dist((esc))))
dev.off()

esc.cor <- cor(t(esc))
epi.cor <- cor(t(epi))

esc.cols <- c(rep("red",29),rep("black",29),rep("blue",29))
epi.cols <- c(rep("red",13),rep("black",13),rep("blue",13))

library(gplots)
png("epi_cor_heatmap.png",width=2000,height=2000)
heatmap.2(epi.cor,scale="none",col=bluered(100),trace="none",mar=c(30,2),labRow="",cexCol=1.8,ColSideColors=epi.cols)
dev.off()

#Remove everything but TF and H3K27ac from esc table
tomatch <- c("H3K27me", "5HMC", "EZH2", "FAIRE", "GFP","5FC","Input")
exclude <- as.list(1:length(tomatch))
for(i in 1:length(tomatch)){
exclude[[i]] <- grep(tomatch[i],rownames(epi))
}
exclude <- sort(unlist(exclude))
rownamestoexclude <- rownames(epi)[exclude]
#esc.tf <- esc[which(!rownames(esc) %in% rownamestoexclude),]
#esc.cols.tf <- esc.cols[which(!rownames(esc) %in% rownamestoexclude)]
#esc.tf.cor <- cor(t(esc.tf))
epi.tf <- epi[which(!rownames(epi) %in% rownamestoexclude),]
epi.cols.tf <- epi.cols[which(!rownames(epi) %in% rownamestoexclude)]
epi.tf.cor <- cor(t(epi.tf))
png("epi_tfonly_cor_heatmap.png",width=2000,height=2000)
heatmap.2(epi.tf.cor,scale="none",col=bluered(100),trace="none",mar=c(30,2),labRow="",cexCol=1.8,ColSideColors=epi.cols.tf)
dev.off()

#Normalise - did not use
library(preprocessCore)
esc.tf.norm <- normalize.quantiles(esc.tf.norm)
rownames(esc.tf.norm) <- rownames(esc.tf)
esc.tf.log10 <- log10(esc.tf+1)

png("esc_tfonly_log10_heatmap.png",width=2000,height=2000)
heatmap.2(esc.tf.log10,scale="none",trace="none",labCol="",mar=c(2,30),cexRow=1.8)
dev.off()
#

#Calculate the mean coverage across subregion and plot 3 column heatmap
esc.ave <- array(NA,dim=c(29,3))
rownames(esc.ave) <- gsub("unmeth.","",rownames(esc)[1:29])
colnames(esc.ave) <- c("Unmeth","DMR","Meth")
for(i in 1:29){
esc.ave[i,1] <- mean(esc[i,])
}
for(i in 1:29){
tmp <- esc[30:58,]
esc.ave[i,2] <- mean(tmp[i,])
}
for(i in 1:29){
tmp <- esc[59:87,]
esc.ave[i,3] <- mean(tmp[i,])
}
esc.tf.ave <- esc.ave[which(!rownames(esc.ave) %in% gsub("unmeth.","",rownamestoexclude[1:9])),]
png("esc_tf_unmethvsmeth_woH3K27ac_heatmap.png",width=2000,height=2000)
heatmap.2(esc.tf.ave3,scale="row",trace="none",col=bluered(100),mar=c(10,5,5,42),cexRow=3,cexCol=3)
dev.off()

epi.ave <- array(NA,dim=c(13,3))
rownames(epi.ave) <- gsub("unmeth.","",rownames(epi)[1:13])
colnames(epi.ave) <- c("Unmeth","DMR","Meth")
for(i in 1:13){
epi.ave[i,1] <- mean(epi[i,])
}
for(i in 1:13){
tmp <- epi[14:26,]
epi.ave[i,2] <- mean(tmp[i,])
}
for(i in 1:13){
tmp <- epi[27:39,]
epi.ave[i,3] <- mean(tmp[i,])
}
epi.tf.ave <- epi.ave[which(!rownames(epi.ave) %in% gsub("unmeth.","",rownamestoexclude[1:5])),]
png("epi_tf_MUvsDMRvsPM_heatmap.png",width=2000,height=2000)
heatmap.2(epi.tf.ave,scale="row",trace="none",col=bluered(100),mar=c(10,50),cexRow=3,cexCol=3)
dev.off()
