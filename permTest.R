# simple permutation test for average coverage within set of regions, compared with random set of regions of same number and size

test.bed <- read.table("ESC SE locations.bed",sep="\t",head=F)
setwd("/data/emmabell42/Coverage_plots/BAM_files")
esrrb <- removeSRduplicates("GSE11431_Esrrb.bam")
ezh2 <- removeSRduplicates("GSE23943_mESC-serum_EZH2.bam")

#first, compute 'statistic'

#for(i in 1:nrow(inputs)){
bam <- get(inputs[i,1])
readcounts.byRegion <- apply(test.bed,1,function(x)sum(chromosome(bam)==as.character(x)[1] & position(bam)<as.numeric(as.character(x)[3]) & (position(bam)+width(bam))<as.numeric(as.character(x)[2])))
#}

# alternative using GRanges

test.bam <- ezh2

reads.GR <- GRanges(seqnames=chromosome(test.bam),ranges=IRanges(start=position(test.bam),end=position(test.bam)+width(test.bam)))
testregions.GR <- GRanges(seqnames=as.character(test.bed[[1]]),ranges=IRanges(start=test.bed[[2]],end=test.bed[[3]]))
overlappingreads <- findOverlaps(subject=testregions.GR,query=reads.GR)
readcounts.total <- length(overlappingreads@queryHits)
readcount.av <- readcounts.total/length(testregions.GR)

# generate random set of regions with same properties as 'test.bed'

all.chrs <- as.character(unique(seqnames(reads.GR)))
chr.max <- sapply(all.chrs,function(x)max(end(reads.GR[seqnames(reads.GR)==x])))

nPerm <- 100
perm.counts <- rep(NA,nPerm)
for(i in 1:nPerm){
# get chromosome lengths
# sample across chromosome lengths
randomGR.maxlen <- sapply(seqnames(testregions.GR),function(x)chr.max[all.chrs==x])
randomGR.maxlen <- randomGR.maxlen-width(testregions.GR)
randomGR.starts <- sapply(randomGR.maxlen,function(x)sample(x,1))
randomGR <- GRanges(seqnames=seqnames(testregions.GR),ranges=IRanges(start=randomGR.starts,end=randomGR.starts+width(testregions.GR)))

random.overlaps <- findOverlaps(subject=randomGR,query=reads.GR)
randomreadcounts.total <- length(random.overlaps@queryHits)
randomreadcount.av <- randomreadcounts.total/length(randomGR)
perm.counts[i] <- randomreadcount.av
}

test.pval <- (sum(perm.counts>readcount.av)+1)/(nPerm+1)
