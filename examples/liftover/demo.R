library(GenomicRanges)

print(load("levine_methyl.RData"))   # anno.levine2
dim(anno.levine2)   # 513 8

anno.levine2$bed <- with(anno.levine2, sprintf("chr%s:%d-%d", CHR, Coordinate_36, Coordinate_36))
tbl <- anno.levine2  # for convenience

write(tbl$bed, file="anno.levine2.bedStrings.txt")

# https://genome.ucsc.edu/cgi-bin/hgLiftOver
# choose original assembly (hg19) and new assembly (hg38)
# Choose File then Submit File
# 504 successful conversion, 9 failed.
# Deleted in new:  chr15:21441686-21441686
# Deleted in new:  chr7:100595769-100595769
# Deleted in new:  chr17:24393906-24393906
# Deleted in new:  chr19:62912474-62912474
# Deleted in new:  chr3:198078873-198078873
# Deleted in new:  chr22:15981381-15981381
# Deleted in new:  chr18:17574536-17574536
# Deleted in new:  chr17:24302936-24302936
# Deleted in new:  chr16:43024-43024
# copy the link address (right-click on "View Conversions"
# curl -O https://genome.ucsc.edu/trash/hglft_genome_70d6_b7db0.bed

txt.hg38 <- scan("hglft_genome_70d6_b7db0.bed", what=character(0), sep="\n")
failures <- c("chr15:21441686-21441686",
              "chr7:100595769-100595769",
              "chr17:24393906-24393906",
              "chr19:62912474-62912474",
              "chr3:198078873-198078873",
              "chr22:15981381-15981381",
              "chr18:17574536-17574536",
              "chr17:24302936-24302936",
              "chr16:43024-43024")



deleters <- match(failures, tbl$bed)
if(length(deleters) > 0)
   tbl <- tbl[-deleters,]

tbl$hg38 <- txt.hg38
tokens <- strsplit(tbl$hg38, "-")
loc.hg38 <- unlist(lapply(tokens, "[", 2))
stopifnot(length(loc.hg38) == nrow(tbl))
tbl$hg38 <- loc.hg38

# prepare for GenomicRanges find overlaps
tbl$chrom <- paste("chr", tbl$CHR, sep="")
tbl$start <- tbl$hg38
tbl$end <- tbl$hg38
preferred.column.order <- c("chrom", "start", "end", "IlmnID", "Name", "Genome_Build", "CHR",
                            "Coordinate_36", "UCSC_RefGene_Name", "Phantom", "Enhancer", "bed", "hg38")
tbl <- tbl[, preferred.column.order]
tbl$start <- as.integer(tbl$start)
tbl$end <- as.integer(tbl$end)

print(load("~/github/genomescale_scripts/geneHancer/v4.7/tbl.enhancers.allGenes.RData"))
head(tbl.enhancers)
head(tbl)
gr.enhancers <- with(tbl.enhancers, GRanges(seqnames=chrom, IRanges(start=start, end=end)))
gr.methylation <- with(tbl, GRanges(seqnames=chrom, IRanges(start=start, end=end)))
tbl.ov <- as.data.frame(findOverlaps(gr.methylation, gr.enhancers))
dim(tbl.ov)
colnames(tbl.ov) <- c("methylationMarks", "enhancers")
tbl.methInEnhancers <- cbind(tbl[tbl.ov$methylationMarks,], tbl.enhancers[tbl.ov$enhancers,])
tbl.methInEnhancers <- tbl.methInEnhancers[order(tbl.methInEnhancers$combinedScore, decreasing=TRUE),]
save(tbl.methInEnhancers, file="tbl.methylationMarksInEnhancers.RData")
