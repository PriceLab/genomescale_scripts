if(!exists("tbl.geneScoresAll")){
   f <- "enhancer_gene_scores.txt"
   tbl.geneScoresAll <- read.table(f, sep="\t", as.is=TRUE, header=TRUE)  #  934287      9
   }
if(!exists("tbl.locsAll")){
   f <- "enhancer_elite_ids.txt"
   tbl.locsAll <- read.table(f, sep="\t", as.is=TRUE, header=TRUE)        #  243281      6
   }

getAllEnhancers <- function(geneSymbol) {
  tbl.geneScores <- subset(tbl.geneScoresAll, symbol==geneSymbol)
  if(nrow(tbl.geneScores) == 0) return(data.frame())
  clusterIDs <- sort(unique(tbl.geneScores$cluster_id))
  tbl.locs <- subset(tbl.locsAll, cluster_id %in% clusterIDs)
  tbl.regions <- merge(tbl.locs, tbl.geneScores, by="cluster_id")
  if(!grepl("chr", tbl.regions$chr[1]))
      tbl.regions$chr <- paste("chr", tbl.regions$chr, sep="")

  table(tbl.regions$regulatory_element_type) # 4 promoters, 21 enhancers

  coi <- c("chr", "enhancer_start", "enhancer_end", "regulatory_element_type", "combined_score")
  tbl.promoters <- subset(tbl.regions, regulatory_element_type=="Promoter/Enhancer")[, coi]
  tbl.enhancers <- subset(tbl.regions, regulatory_element_type=="Enhancer")[, coi]

  #track.p <- DataFrameAnnotationTrack("GH 4.7 promoters", tbl.promoters, color="blue", trackHeight=30)
  #track.e <- DataFrameAnnotationTrack("GH 4.7 enhancers", tbl.enhancers, color="darkblue", trackHeight=30)
  #displayTrack(igv, track.p)
  #displayTrack(igv, track.e)
  #region("enhancers")

  # list(promoters=tbl.promoters, enhancers=tbl.enhancers)
  tbl.out <- tbl.regions[, coi]
  tbl.out$geneSymbol <- geneSymbol
  tbl.out

} # getAllEnhancers
