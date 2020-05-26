library(TrenaProjectHG38.generic)   # for access to hg38 genehancer
library(GenomicRanges)
library(rtracklayer)

#------------------------------------------------------------
# get all TREM2 enhancers, then get that slice of the bigwig
#------------------------------------------------------------
tp <- TrenaProjectHG38.generic()
tbl.gh <- getEnhancers(tp, "TREM2")

gr.region <- with(tbl.gh, GRanges(seqnames=chrom[1], IRanges(start=min(start)-5000,
                                                             end=max(end)+5000)))
#------------------------------------------------------------
# my randomly selected bigwig file is here
#   scp pshannon@khaleesi:s/data/public/GSM2579603/DLPFC_neuron.bw .
#------------------------------------------------------------

your.bigwig.file <- "DLPFC_neuron.bw"
gr.small <- import(con=your.bigwig.file, which=gr.region)
length(gr.small)
#------------------------------------------------------------
# keep just the high-scoring peaks
#------------------------------------------------------------
fivenum(gr.small$score)  #  1  1  2  3 19
gr.strongPeaks <- gr.small[gr.small$score >= 7] # 60

#----------------------------------------------------------------------
# find the overlaps.  GRanges knows how to convert the gh data.frame
#----------------------------------------------------------------------

tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.gh), gr.strongPeaks, type="any"))
colnames(tbl.ov) <- c("gh", "bw")
dim(tbl.ov)    # 14 2

length(unique(tbl.ov$gh)) # 4 gh regions
length(unique(tbl.ov$bw)) # 14 peaks

#----------------------------------------------------------------------
# combine the two tables, column-wise.
#----------------------------------------------------------------------
tbl.peaksInEnhancers <- cbind(tbl.gh[tbl.ov$gh,], gr.strongPeaks[tbl.ov$bw])
dim(tbl.peaksInEnhancers)  # 14 22


