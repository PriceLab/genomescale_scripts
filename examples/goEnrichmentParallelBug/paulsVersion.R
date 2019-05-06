library(GOstats)
library(GO.db)
library(Category)
library(org.Hs.eg.db)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl.models")){
   tbl.models <- get(load("lympho-footprints-2019feb22-tbl.models.all.RData"))
   dim(tbl.models)
   printf("tfs with models: %d", length(unique(tbl.models$tf.symbol)))
   }

#---------------------------------------------------------------------------------
source("https://raw.githubusercontent.com/PriceLab/genomescale_scripts/master/symToGeneID.R"); test_assignGeneIDs()
#------------------------------------------------------------------------------------------------------------------------
goEnrichment <- function(geneSymbols)
{
  symbol.entrez.map <- assignGeneIDs(geneSymbols)
  gene.universe = character(0)
  geneIDs <- unlist(symbol.entrez.map$mapped, use.names=FALSE)

  go.params <- new("GOHyperGParams", geneIds=unique(geneIDs),
                   universeGeneIds = gene.universe, annotation = "org.Hs.eg.db",
                   ontology = 'BP', pvalueCutoff = 0.05, conditional = FALSE,
                   testDirection = "over")

  go.bp.hgr <- hyperGTest(go.params)
  tbl.go <- summary(go.bp.hgr)
  geneSymbols <- lapply(tbl.go$GOBPID,
                        function(goTerm){
                          all.geneIDs.this.term <- unique(unlist(get(goTerm, org.Hs.egGO2ALLEGS)))
                          keepers <- intersect(geneIDs, all.geneIDs.this.term)
                          keeper.geneSymbols <- mget(keepers, org.Hs.egSYMBOL)
                          keeper.geneSymbols <- unlist(keeper.geneSymbols, use.names=FALSE)
                          paste(keeper.geneSymbols, collapse=";")
                        })
  #RSQLite::dbDisconnect(dbconn(go.params@datPkg@db))
  tbl.go$genes <- unlist(geneSymbols, use.names=FALSE)
  rm(go.params)
  tbl.go

} # goEnrichment
#------------------------------------------------------------------------------------------------------------------------
test_goEnrichment <- function(junk)
{
   printf("--- test_goEnrichment")

   goi <- subset(tbl.models, tf.symbol == "TBX15" & rank <= 2 & abs(spearmanCoeff) > 0.5)$target.symbol
   length(goi)
   tbl.go <- goEnrichment(goi)
   checkTrue(is.data.frame(tbl.go))
   checkEquals(ncol(tbl.go), 8)
   checkTrue(nrow(subset(tbl.go, Pvalue < 10e-5)) > 0)

} # test_goEnrichment
#------------------------------------------------------------------------------------------------------------------------

