library(MotifDb)
library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP150.GRCh38)

snps <- c("rs483082", "rs59325138", "rs584007", "rs438811")
demo.snp <- snps[4]

apoc1.tfs <- c("SPI1", "IRF5", "ELK3", "STAT3", "IKZF1", "RREB1", "JUNB", "SMAD2", "FOS", "FLI1")

matches <- unlist(lapply(apoc1.tfs,
                         function(tf) grep(sprintf("^%s$", tf),mcols(MotifDb)$geneSymbol)))

mini.db <- MotifDb[matches]
motifs <- query(mini.db, c("hsapiens", "jaspar2018"))
length(motifs)

snps.gr <- snps.from.rsid(rsid = demo.snp,
                          dbSNP=SNPlocs.Hsapiens.dbSNP150.GRCh38,
                          search.genome=BSgenome.Hsapiens.UCSC.hg38)

results <- motifbreakR(snpList = snps.gr,
                       filterp = TRUE,
                       pwmList = motifs,
                       show.neutral=FALSE,
                       method = c("ic", "log", "notrans")[2],
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::bpparam(),

                       verbose=TRUE)
length(results)  # 5
length(results[results$effect=="strong"])  # 4

plotMB(results, rsid=demo.snp, effect=c("strong")) # takes a minute or two



results <- calculatePvalue(results)

  # translate the results into a more convenient data structure, a data.frame
  # add a few summary columns.

results$rsid <- names(results)
tbl.breaks <- as.data.frame(results, row.names=seq_len(length(results)))
tbl.breaks$delta <- with(tbl.breaks, pctRef-pctAlt)
tbl.breaks$refPvalScore <- -log10(tbl.breaks$Refpvalue)
tbl.breaks$altPvalScore <- -log10(tbl.breaks$Altpvalue)
tbl.breaks$deltaPvalScore <- tbl.breaks$refPvalScore - tbl.breaks$altPvalScore

  # plot the results




#                                rs438811
# seqnames                          chr19
# start                          44913479
# end                            44913488
# width                                10
# strand                                -
# REF                                   C
# ALT                                   T
# snpPos                         44913484
# motifPos                              5
# geneSymbol                  FOSL2::JUNB
# dataSource                   jaspar2018
# providerName                   MA1138.1
# providerId                     MA1138.1
# seqMatch     ttgaaCtttt
# pctRef                        0.4863753
# pctAlt                        0.6256244
# scoreRef                      -19.41419
# scoreAlt                      -11.11889
# Refpvalue                          <NA>
# Altpvalue                          <NA>
# alleleRef                             0
# alleleAlt                             1
# effect                           strong


#head(results)
# GRanges object with 6 ranges and 18 metadata columns:
#              seqnames            ranges strand |            REF            ALT    snpPos  motifPos  geneSymbol  dataSource providerName  providerId                            seqMatch            pctRef            pctAlt          scoreRef          scoreAlt Refpvalue Altpvalue         alleleRef         alleleAlt      effect
#                 <Rle>         <IRanges>  <Rle> | <DNAStringSet> <DNAStringSet> <integer> <integer> <character> <character>  <character> <character>                         <character>         <numeric>         <numeric>         <numeric>         <numeric> <logical> <logical>         <numeric>         <numeric> <character>
#     rs438811    chr19 44913479-44913488      + |              C              T  44913484         6        ALX3  jaspar2018     MA0634.1    MA0634.1          ttgaaCtttt                 0.601057413561636 0.768198052825913 -2.57027506727742   2.1219726731512      <NA>      <NA> 0.008963027511515 0.980580107058384      strong
#     rs438811    chr19 44913470-44913484      + |              C              T  44913484        15          AR  jaspar2018     MA0007.2    MA0007.2 ggctaatttttgaaC                     0.688128672731351 0.674006861629654 -11.5470109572363 -12.6745490941078      <NA>      <NA> 0.627520970908442 0.203194717115831        weak
#     rs438811    chr19 44913481-44913488      - |              C              T  44913484         5 ARNT::HIF1A  jaspar2018     MA0259.1    MA0259.1            gaaCtttt                 0.614729466124871 0.446963017109685 -5.69137723905057 -11.7244634608494      <NA>      <NA>                 1                 0      strong
#     rs584007    chr19 44913217-44913224      + |              A              G  44913221         5 ARNT::HIF1A  jaspar2018     MA0259.1    MA0259.1           ataaAagc                  0.477961455515845 0.645727904531031 -10.6097217902514 -4.57663556845258      <NA>      <NA>                 0                 1      strong
#   rs59325138    chr19 44913029-44913036      + |              C              T  44913034         6 ARNT::HIF1A  jaspar2018     MA0259.1    MA0259.1          ggttgCgg                   0.568634715349585 0.736401164364772 -7.34900072355583 -1.31591450175703      <NA>      <NA>                 0                 1      strong
#     rs483082    chr19 44912916-44912923      + |              G              T  44912921         6 ARNT::HIF1A  jaspar2018     MA0259.1    MA0259.1          gggctGac                     0.4629530857331 0.630719534748286 -11.1494410814414 -5.11635485964258      <NA>      <NA>                 0                 1      strong


# reached end of SNPs list length = 4 with 1265 potentially disruptive
# matches to 450 of 537 motifs.

# loc <- 88884592
# snp.name <- sprintf("%s:%d:%s:%s", "chr5", loc, "A", "C")
# tbl.snp <- data.frame(chrom="chr5", start=loc-1, end=loc, name=snp.name, score=0, strand="+",
#                       stringsAsFactors=FALSE)
#
# bedFile <- "snp.bed"
# write.table(tbl.snp, file=bedFile, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
# getSeq(BSgenome.Hsapiens.UCSC.hg38, "chr5", loc, loc)
# snps.gr <- snps.from.file(bedFile,
#                           search.genome = BSgenome.Hsapiens.UCSC.hg38,
#                           format = "bed")
#
# motifs <- query(MotifDb, c("sapiens", "TBR1"))
# motifs <- query(MotifDb, c("sapiens", "jaspar2018"), c("EGR3", "TBR1"))
# results <- motifbreakR(snpList = snps.gr,
#                           filterp = TRUE,
#                           pwmList = motifs,
#                           show.neutral=FALSE,
#                           method = c("ic", "log", "notrans")[2],
#                           bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
#                           BPPARAM = BiocParallel::bpparam(),
#                           verbose=TRUE)
# plotMB(results, rsid=snp.name, effect=c("weak", "strong")) # takes a minute or two


