#library(GOstats)
#library(Category)
library(BiocParallel)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tbl.models")){
   tbl.models <- get(load("lympho-footprints-2019feb22-tbl.models.all.RData"))
   dim(tbl.models)
   printf("tfs with models: %d", length(unique(tbl.models$tf.symbol)))
   }

#---------------------------------------------------------------------------------
source("https://raw.githubusercontent.com/PriceLab/genomescale_scripts/master/symToGeneID.R"); test_assignGeneIDs()
#------------------------------------------------------------------------------------------------------------------------
#goEnrichment <- function(geneSymbols)
goEnrichment <- function(geneIDs)
{
  gene.universe = character(0)

  require("org.Hs.eg.db")
  require("GO.db")
  require("GOstats")
  go.params <- new("GOHyperGParams",
                   geneIds=unique(geneIDs),
                   universeGeneIds = gene.universe,
                   annotation = "org.Hs.eg.db",
                   ontology = 'BP',
                   pvalueCutoff = 0.05,
                   conditional = FALSE,
                   testDirection = "over")

  go.bp.hgr <- hyperGTest(go.params)
  tbl.go <- summary(go.bp.hgr)
  #geneSymbols <- lapply(tbl.go$GOBPID,
  #                      function(goTerm){
  #                        all.geneIDs.this.term <- unique(unlist(get(goTerm, org.Hs.egGO2ALLEGS)))
  #                        keepers <- intersect(geneIDs, all.geneIDs.this.term)
  #                        keeper.geneSymbols <- mget(keepers, org.Hs.egSYMBOL)
  #                        keeper.geneSymbols <- unlist(keeper.geneSymbols, use.names=FALSE)
  #                        paste(keeper.geneSymbols, collapse=";")
  #                      })
  RSQLite::dbDisconnect(dbconn(go.params@datPkg@db$conn))
  #tbl.go$genes <- unlist(geneSymbols, use.names=FALSE)
  # browser()
  unloadNamespace("org.Hs.eg.db")
  unloadNamespace("GOstats")
  unloadNamespace("GO.db")
  rm(go.params)
  gc()
  tbl.go

} # goEnrichment
#------------------------------------------------------------------------------------------------------------------------
test_goEnrichment <- function(junk)
{
   printf("--- test_goEnrichment")

   goi <- subset(tbl.models, tf.symbol == "TBX15" & rank <= 2 & abs(spearmanCoeff) > 0.5)$target.symbol
   length(goi)
   symbol.entrez.map <- assignGeneIDs(goi)
   geneIDs <- unlist(symbol.entrez.map$mapped, use.names=FALSE)

   tbl.go <- goEnrichment(goi)
   checkTrue(is.data.frame(tbl.go))
   checkEquals(ncol(tbl.go), 8)
   checkTrue(nrow(subset(tbl.go, Pvalue < 10e-5)) > 0)

} # test_goEnrichment
#------------------------------------------------------------------------------------------------------------------------
if(1 == 2){
   param=MulticoreParam(workers=12, stop.on.error=FALSE)
   #param=MulticoreParam(stop.on.error=FALSE)
   #param=SerialParam()
   max <- 10
   load("test.geneIDs.RData")
   system.time(results <- bplapply(1:max,
                                   function(geneIDs){goEnrichment(geneIDs)},
                                   BPPARAM=param))
   }
   
               

