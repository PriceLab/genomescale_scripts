library(TrenaProjectHG38.generic)   # for access to hg38 genehancer
library(GenomicRanges)
library(rtracklayer)
library(igvR)
library(RPostgreSQL)
#----------------------------------------------------------------------------------------------------
if(!exists("igv")){
   igv <- igvR()
   setGenome(igv, "hg38")
   }
#----------------------------------------------------------------------------------------------------
geneHancerByGene <- function(goi)
{
  showGenomicRegion(igv, goi)
  tp <- TrenaProjectHG38.generic()
  tbl.gh <- getEnhancers(tp, goi, maxSize=100000)
  loc <- getGenomicRegion(igv)

  shoulder <- 1000000
  loc$start <- loc$start - shoulder
  loc$end <- loc$end + shoulder

  with(tbl.gh, showGenomicRegion(igv, sprintf("%s:%d-%d", chrom[1],
                                              min(start)-shoulder,
                                             max(end)+shoulder)))

  track.name <- sprintf("%s-GH", goi)
  track <- DataFrameQuantitativeTrack("GeneHancer", tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                                      autoscale=TRUE, color="brown")
  displayTrack(igv, track)

}
#----------------------------------------------------------------------------------------------------
initializeIgv <- function(gene)
{
   showGenomicRegion(igv, gene)
   loc <- getGenomicRegion(igv)

   shoulder <- 1000000
   loc$start <- loc$start - shoulder
   loc$end <- loc$end + shoulder
   with(loc, showGenomicRegion(igv, sprintf("%s:%d-%d", chrom, start, end)))

} # initializeIgv
#----------------------------------------------------------------------------------------------------
createAndDisplayTrack <- function(bw.file, trackName)
{
   stopifnot(file.exists(bw.file))
   roi <- getGenomicRegion(igv)
   groi <- with(roi, GRanges(seqnames=chrom, IRanges(start=start, end=end)))
   gr.small <- import(con=bw.file, which=groi)
   length(gr.small) # 25088
   tbl.chip <- as.data.frame(gr.small)[, c("seqnames", "start", "end", "score")]
   colnames(tbl.chip)[1] <- "chrom"
   tbl.chip$chrom <- as.character(tbl.chip$chrom)
   track <- DataFrameQuantitativeTrack(trackName, tbl.chip, autoscale=TRUE, color="random")
   igvR::displayTrack(igv, track)
   invisible(tbl.chip)

} # createAndDisplayTrack
#----------------------------------------------------------------------------------------------------
getGeneHancerByRegion <- function(chrom, start, end)
{
    db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="gh411", host="khaleesi")


    query <- paste0("select e.chr as chrom, ",
                    "e.element_start as start, ",
                    "e.element_end as end, ",
                    "a.symbol as gene, ",
                    "a.eqtl_score as eqtl, ",
                    "a.chic_score as HiC, ",
                    "a.erna_score as erna, ",
                    "a.expression_score as coexpression, ",
                    "a.distance_score as distanceScore, ",
                    "a.tss_proximity as tssProximity, ",
                    "a.combined_score as combinedScore, ",
                    "a.is_elite as elite, ",
                    "t.source as source, ",
                    "t.tissue as tissue, ",
                    "e.type as type, ",
                    "a.ghid as ghid ",
                    "from associations AS a, ",
                    "tissues AS t, elements as e ",
                    "where e.chr='%s' and e.element_start > %d and e.element_end < %d",
                    "AND a.ghid=t.ghid ",
                    "AND e.ghid=a.ghid")
    query <- sprintf(query, chrom, start, end)
    tbl <- dbGetQuery(db, query)
    printf("nrow tbl.gh: %d", nrow(tbl))
    #tbl <- subset(tbl, combinedscore > 5)
    dim(tbl)
    dbDisconnect(db)

    invisible(tbl)

} # getGeneHancerByRegion
#----------------------------------------------------------------------------------------------------
displayChip <- function()
{
   data.directory <- "~/github/genomescale_scripts/bobRostomily/AM.HMRI.TWIST1.CSeq.37738.2"
   bw.files <- list(
      twE12="1_079L_00YMHMRI_T98G-TW-E12_TWIST1_hg38_i88_dmnorm_signal.bw",
      twS68AE="2_079M_00YMHMRI_T98G-TW-S68AE_TWIST1_hg38_i89_dmnorm_signal.bw",
      twtw="3_079N_00YMHMRI_T98G-TW-TW_TWIST1_hg38_i90_dmnorm_signal.bw",
      twTwist1="4_079O_00YMHMRI_T98G-TW_TWIST1_hg38_i91_dmnorm_signal.bw")
      #pooled="5_079T_00YMHMRI_Pooled_Input_hg38_i92_uniqnorm_signal.bw")

    tbls.chip <- list()

    bw.files <- bw.files[1]

    for(f in seq_len(length(bw.files))){
       trackName <- names(bw.files)[f]
       file <- file.path(data.directory, bw.files[[f]])
       tbl.chip <- createAndDisplayTrack(file, trackName)
       tbls.chip[[f]] <- subset(tbl.chip, score > 10)
       } # for f

    names(tbls.chip) <- names(bw.files)
    invisible(tbls.chip)

} # displayChip
#----------------------------------------------------------------------------------------------------
intersectBigChipWithRegionalGeneHancer <- function(bigChip)
{
   roi <- getGenomicRegion(igv)
   tbl.gh <- with(roi, getGeneHancerByRegion(chrom, start, end))
   tbl.gh <- subset(tbl.gh, combinedscore >= 10)
   tbl.gh <- subset(tbl.gh, nchar(gene) < 10)    # try to get just the gene symbols
   dim(tbl.gh)
   unique(tbl.gh$gene)

   tbl.ov.1 <- as.data.frame(findOverlaps(GRanges(tbl.gh), GRanges(subset(bigChip[[1]], score > 10))))
   colnames(tbl.ov.1) <- c("gh", "chip")
   tbl.gh.chip <- tbl.gh[unique(tbl.ov.1$gh),]

   for(targetGene in unique(tbl.gh.chip$gene)){
      tbl.sub <- subset(tbl.gh.chip, gene==targetGene)
      track.name <- sprintf("%s-%s-%s", names(bigChip)[1], "gh", targetGene)
      track <- DataFrameAnnotationTrack(track.name,
                                        unique(tbl.sub[, c("chrom", "start", "end", "gene")]),
                                        trackHeight=22,
                                        color="random")
      displayTrack(igv, track)
      } # for targetGene

   track.name <- sprintf("%s-%s", names(bigChip)[1], "ChIP")
   track <- DataFrameAnnotationTrack(track.name,
                                     unique(bigChip[[1]][, c("chrom", "start", "end")]),
                                     color="darkred")
   displayTrack(igv, track)

   findOverlaps(GRanges(tbl.gh), GRanges(bigChip[[2]]))
   findOverlaps(GRanges(tbl.gh), GRanges(bigChip[[3]]))
   findOverlaps(GRanges(tbl.gh), GRanges(bigChip[[4]]))
   dim(tbl.ov)

} # intersectBigChipWithRegionalGeneHancer
#----------------------------------------------------------------------------------------------------
run <- function()
{
   goi <- "PDGFRA"
   geneHancerByGene(goi)
   bigChip <- displayChip()

} # run
#----------------------------------------------------------------------------------------------------
