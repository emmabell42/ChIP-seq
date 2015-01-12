library(GenomicRanges)
gal <- readGAlignments("/path/to/your/file.bam")
as(gal, "GRanges")
