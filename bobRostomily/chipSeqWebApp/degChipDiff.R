library(EnsDb.Hsapiens.v79)
library(TrenaProjectHG38.generic)   # for access to hg38 genehancer
library(GenomicRanges)

deg.directory <-   "~/github/genomescale_scripts/bobRostomily/chipSeqWebApp/deg"
deg.file <- file.path(deg.directory, "twist_E12_dimer-vs-WT.tsv")
file.exists(deg.file)
tbl.deg <- read.table(deg.file, header=TRUE, as.is=TRUE, quote=NULL, nrow=-1)

fivenum(tbl.deg$log2FoldChange)   # -7.368000 -0.509415  0.190860  0.720675 10.651000
fivenum(-log10(tbl.deg$pvalue))   # 1.702677  3.217786  6.848416 19.072194       Inf

# find just a few genes to experiment with

ensg.deg.pos <- subset(tbl.deg, log2FoldChange >  8 & -log10(pvalue) > 8)$Gene_ID
ensg.deg.neg <- subset(tbl.deg, log2FoldChange < -5 & -log10(pvalue) > 8)$Gene_ID

sym.deg.pos <- select(EnsDb.Hsapiens.v79, key=ensg.deg.pos,  columns=c("SYMBOL", "GENEID"),  keytype="GENEID")$SYMBOL
sym.deg.neg <- select(EnsDb.Hsapiens.v79, key=ensg.deg.neg,  columns=c("SYMBOL", "GENEID"),  keytype="GENEID")$SYMBOL

tp <- TrenaProjectHG38.generic()
tbl.peaks.twe12 <- read.table("twe12.bed", sep="\t", skip=1, header=FALSE, as.is=TRUE)[, c(1,2,3,5)]
tbl.peaks.twtw <- read.table("twtw.bed", sep="\t", skip=1, header=FALSE, as.is=TRUE)[, c(1,2,3,5)]
colnames(tbl.peaks.twe12) <- c("chrom", "start", "end", "score")
colnames(tbl.peaks.twtw) <- c("chrom", "start", "end", "score")
gr.peaks.twe12 <- GRanges(tbl.peaks.twe12)
gr.peaks.twtw <- GRanges(tbl.peaks.twtw)

gene <- sym.deg.pos[1]
gene <- "POSTN"

igv <- start.igv(gene)

for(gene in sym.deg.pos){
   setBrowserWindowTitle(igv, gene)
   showGenomicRegion(igv, gene)
   tbl.gh <- getEnhancers(tp, gene, maxSize=40000)
   tbl.gh <- subset(tbl.gh, elite)
   shoulder <- 1000
   with(tbl.gh, showGenomicRegion(igv,
                         sprintf("%s:%d-%d", chrom[1], min(start-shoulder), max(end) + shoulder)))
   track <- DataFrameQuantitativeTrack("gh",
                                       tbl.gh[, c("chrom", "start", "end", "combinedscore", "gene", "ghid")],
                                       autoscale=FALSE, min=0, max=50, color="brown")
   displayTrack(igv, track)
   gr.gh <- GRanges(tbl.gh)
   tbl.ov.twe12 <- as.data.frame(findOverlaps(gr.peaks.twe12, gr.gh,type="any"))
   tbl.ov.twtw <- as.data.frame(findOverlaps(gr.peaks.twtw,   gr.gh,type="any"))

   tbl.twe12.hits <- as.data.frame(gr.peaks.twe12[tbl.ov.twe12$queryHits])[, c(1,2,3,6)]
   tbl.twe12.hits$seqnames <-  as.character(tbl.twe12.hits$seqnames)
   track <- DataFrameQuantitativeTrack("twe12", tbl.twe12.hits, autoscale=TRUE,color="darkgreen")
   displayTrack(igv, track)

   tbl.twtw.hits <- as.data.frame(gr.peaks.twtw[tbl.ov.twtw$queryHits])[, c(1,2,3,6)]
   tbl.twtw.hits$seqnames <-  as.character(tbl.twtw.hits$seqnames)
   track <- DataFrameQuantitativeTrack("twtw", tbl.twtw.hits, autoscale=TRUE,color="red")
   displayTrack(igv, track)

   invisible(readline(prompt="Press [enter] to continue: "))
   }




