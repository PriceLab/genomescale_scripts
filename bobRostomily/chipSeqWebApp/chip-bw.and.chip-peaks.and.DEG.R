# look for correlations between gene expression,
# twist_deletion-vs-WT.tsv
#
#             original name                   my preferred name
#  TEvsCNTRL T98G_DEG_all.tsv         ../deg/twist_E12_dimer-vs-WT
#  TS68AEvsCNTRL T98G_DEG_all.tsv     ../deg/ts68ae_mutant-vs-WT
#  TS68AEvsTE T98G DEG_all.tsv        ../deg/ts68ae_mutant-vs-twist_E12_dimer
#  TWmonomer vsCNTRLWT_DEG_all.tsv    ../deg/twistMonomer-vs-T98G
#  TwTwvsCNTRL T98G_DEG_all.tsv       ../deg/twist_twist_dimer-vs-WT
#  dTWB98vsCNTRL98_DEG_all.tsv        ../deg/twist_deletion-vs-WT
#----------------------------------------------------------------------------------------------------
library(EnsDb.Hsapiens.v79)
library(igvR)


deg.directory <- "~/github/genomescale_scripts/bobRostomily/chipSeqWebApp/deg"

chipseq.directory <- "~/github/genomescale_scripts/bobRostomily/AM.HMRI.TWIST1.CSeq.37738.2"

bw.files <- list(
     twE12="1_079L_00YMHMRI_T98G-TW-E12_TWIST1_hg38_i88_dmnorm_signal.bw",
     twS68AE="2_079M_00YMHMRI_T98G-TW-S68AE_TWIST1_hg38_i89_dmnorm_signal.bw",
     twtw="3_079N_00YMHMRI_T98G-TW-TW_TWIST1_hg38_i90_dmnorm_signal.bw",
     twTwist1="4_079O_00YMHMRI_T98G-TW_TWIST1_hg38_i91_dmnorm_signal.bw")

deg.file <- file.path(deg.directory, "twist_E12_dimer-vs-WT.tsv")
file.exists(deg.file)
chip.file.01 <- file.path(chipseq.directory, bw.files[["twE12"]])
chip.file.02 <- file.path(chipseq.directory, bw.files[["twtw"]])

chipDiffsFile <- file.path("~/github/genomescale_scripts/bobRostomily/AM.HMRI.TWIST1.CSeq.37738.2",
                           "00YMHMRI_TWIST1_mergedregs_topdiffs_TW-E12vsTW-TW.txt")
file.exists(chipDiffsFile)
tbl.chipDiffs <- read.table(chipDiffsFile, sep="\t", as.is=TRUE, nrow=-1, quote=NULL, fill=TRUE, header=TRUE)
dim(tbl.chipDiffs)  # 981 33

tbl <- read.table(deg.file, sep="\t", as.is=TRUE, header=TRUE, quote=NULL)
dim(tbl)  # 9263 30
fivenum(tbl$log2FoldChange)  # -7.368000 -0.509415  0.190860  0.720675 10.651000
fivenum(-log10(tbl$pvalue))  # 1.702677  3.217786  6.848416 19.072194       Inf

tbl.nz <- subset(tbl, pvalue > 0)
with(tbl.nz, plot(log2FoldChange, -log10(pvalue)))
tbl.super.deg <- subset(tbl.nz, abs(log2FoldChange) > 5 & -log10(pvalue) > 100)
dim(tbl.super.deg)
ensg.oi <- tbl.super.deg$Gene_ID #  "ENSG00000116774" "ENSG00000173546"

#sym.oi <- select(EnsDb.Hsapiens.v79, key=ensg.oi,  columns=c("SYMBOL", "GENEID"),  keytype="GENEID")$SYMBOL
 # "OLFML3" "CSPG4"

 # "DCN" "MAST4" "RUNX1T1" "NFATC4" "PRUNE2" "MAN1A1" "PHACTR1" "C7"
 # "OLFML3" "PPP1R3C" "CSMD2" "GALNT5" "SULF1" "ADAMTS14" "FGF5" "LUM"
 # "FBLN5" "TINAGL1" "GFRA1" "PRDM8" "ACTG2" "SGCD" "LRRC15" "RCAN2"
 # "CSPG4" "MYO1D" "SCN4B" "GAS1"

sym.oi <- c("POSTN") # , "PDGFRA")
library(TrenaProjectHG38.generic)   # for access to hg38 genehancer
library(GenomicRanges)
library(rtracklayer)
        #------------------------------------------------------------
        # get all TREM2 enhancers, then get that slice of the bigwig
        #------------------------------------------------------------
tp <- TrenaProjectHG38.generic()


for(geneSymbol in sym.oi){
   tbl.gh <- getEnhancers(tp, geneSymbol, maxSize=40000)
   gr.gh <- GRanges(tbl.gh)
   gr.chip.01 <- import(con=chip.file.01, which=gr.gh)
   gr.chip.02 <- import(con=chip.file.02, which=gr.gh)
   tbl.ov.01 <- as.data.frame(findOverlaps(gr.gh, gr.chip.01))
   colnames(tbl.ov.01) <- c("gh", "chip.01")
   tbl.ov.02 <- as.data.frame(findOverlaps(gr.gh, gr.chip.02))
   colnames(tbl.ov.02) <- c("gh", "chip.02")
   printf("-- %s: %d  %d", geneSymbol, length(gr.chip.01), length(gr.chip.02))
   }

igv <- start.igv("POSTN")
track <- DataFrameQuantitativeTrack("gh", tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                                    autoscale=TRUE, color="brown")
displayTrack(igv, track)


gr.region <- with(tbl.gh, GRanges(seqnames=chrom[1], IRanges(start=min(start)-5000,
                                                             end=max(end)+ 5000)))
gr.chip <- import(con=chip.file.01, which=gr.region)
length(gr.chip)
tbl.chip <- as.data.frame(gr.chip)[, c("seqnames", "start", "end", "score")]
colnames(tbl.chip)[1] <-"chrom"
tbl.chip$chrom <- as.character(tbl.chip$chrom)
lapply(tbl.chip, class)

track <- DataFrameQuantitativeTrack("twe12", tbl.chip[, c("chrom", "start", "end", "score")],
                                    autoscale=TRUE, color="red")
displayTrack(igv, track)

gr.chip.01 <- import(con=chip.file.01, which=gr.region)
length(gr.chip.01)
tbl.chip <- as.data.frame(gr.chip.01)[, c("seqnames", "start", "end", "score")]
colnames(tbl.chip)[1] <-"chrom"
tbl.chip$chrom <- as.character(tbl.chip$chrom)
lapply(tbl.chip, class)

track <- DataFrameQuantitativeTrack("twe12", tbl.chip[, c("chrom", "start", "end", "score")],
                                    autoscale=TRUE, color="red")
displayTrack(igv, track)



gr.chip.02 <- import(con=chip.file.02, which=gr.region)
length(gr.chip.02)
tbl.chip.02 <- as.data.frame(gr.chip.02)[, c("seqnames", "start", "end", "score")]
colnames(tbl.chip.02)[1] <-"chrom"
tbl.chip.02$chrom <- as.character(tbl.chip.02$chrom)
lapply(tbl.chip.02, class)

track <- DataFrameQuantitativeTrack("twtw", tbl.chip.02[, c("chrom", "start", "end", "score")],
                                    autoscale=TRUE, color="blue")
displayTrack(igv, track)


#------------------------------------------------------------
# now the macs2 peaks
#------------------------------------------------------------
# chipPeaks.directory <- "./"
# chipPeak.files <- list(
#    twe12="twe12.bed",
#    twtw="twtw.bed")

tbl.peaks.twe12 <- read.table("twe12.bed", sep="\t", skip=1, header=FALSE, as.is=TRUE)[, c(1,2,3,5)]
tbl.peaks.twtw <- read.table("twtw.bed", sep="\t", skip=1, header=FALSE, as.is=TRUE)[, c(1,2,3,5)]
colnames(tbl.peaks.twe12) <- c("chrom", "start", "end", "score")
colnames(tbl.peaks.twtw) <- c("chrom", "start", "end", "score")

x <- getGenomicRegion(igv)

tbl.sub.twe12 <- subset(tbl.peaks.twe12, chrom==x$chrom & start > x$start & end < x$end)
dim(tbl.sub.twe12)
track <- DataFrameQuantitativeTrack("TWE12-peaks", tbl.sub.twe12, autoscale=FALSE, min=0, max=500, color="random")
displayTrack(igv, track)
tbl.sub.twtw <- subset(tbl.peaks.twtw, chrom==x$chrom & start > x$start & end < x$end)
track <- DataFrameQuantitativeTrack("TWTW-peaks", tbl.sub.twtw, autoscale=FALSE, min=0, max=500, color="random")
displayTrack(igv, track)
