library(TrenaProject)     # has GeneHancer data in single table
library(GenomicRanges)

   # get tbl.enhancers
full.path <- system.file(package="TrenaProject", "extdata", "geneHancer.v4.7.allGenes.RData")
stopifnot(file.exists(full.path))
print(load(full.path))
head(tbl.enhancers)

   # some snps for demonstration
tbl.snps.20 <- read.table("tbl.snps.20.tsv", sep="\t", as.is=TRUE, header=TRUE)

  # create GRanges represeentations of the two tables, for fast intersection
gr.snps <- GRanges(tbl.snps.20)
gr.enhancers <- GRanges(tbl.enhancers)

tbl.ov <- as.data.frame(findOverlaps(gr.snps, gr.enhancers))
dim(tbl.ov)
colnames(tbl.ov) <- c("snp", "enhancer")
dim(tbl.ov)
head(tbl.ov)

  # find all the genes with enhancers containing snp # 2
tbl.snps[2,]
tbl.enhancers[subset(tbl.ov, snp==2)$enhancer,]

