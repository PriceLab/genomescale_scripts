# library(GOstats)
library(BiocParallel)
#----------------------------------------------------------------------------------------------------
goEnrichment <- function(geneSymbols)
{


  source("https://raw.githubusercontent.com/PriceLab/genomescale_scripts/master/symToGeneID.R");
  symbol.entrez.map <- assignGeneIDs(geneSymbols)
  geneIDs <- unlist(symbol.entrez.map$mapped, use.names=FALSE)

  gene.universe = character(0)

  require("GO.db")
  require("GOstats")
  #browser()
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

  geneSymbols <- lapply(tbl.go$GOBPID,
                        function(goTerm){
                          all.geneIDs.this.term <- unique(unlist(get(goTerm, org.Hs.egGO2ALLEGS)))
                          keepers <- intersect(geneIDs, all.geneIDs.this.term)
                          keeper.geneSymbols <- mget(keepers, org.Hs.egSYMBOL)
                          keeper.geneSymbols <- unlist(keeper.geneSymbols, use.names=FALSE)
                          paste(keeper.geneSymbols, collapse=";")
                        })
  tbl.go$genes <- unlist(geneSymbols, use.names=FALSE)

  on.exit({RSQLite::dbDisconnect(dbconn(GO.db));
           RSQLite::dbDisconnect(go.params@datPkg@db$conn);
           unloadNamespace("GOstats");
           unloadNamespace("GO.db");
           unloadNamespace("Category");
           unloadNamespace("org.Hs.eg.db");
           rm(go.params);
           gc();
           })
  tbl.go

} # goEnrichment
#----------------------------------------------------------------------------------------------------
run <- function() {
   geneIDs <- get(load("test.geneIDs.RData"))[1:10]
   geneSymbols <- get(load("test.geneSymbols.RData"))[1:10]
   length(geneIDs)
   length(geneSymbols)
   x <- goEnrichment(geneIDs)
   x <- goEnrichment(geneSymbols)
   max <- 40
   param=MulticoreParam(workers=40, stop.on.error=FALSE)
   system.time(results <- bplapply(1:max, function(i){goEnrichment(geneIDs)}, BPPARAM=param))
   system.time(results <- bplapply(1:max, function(i){goEnrichment(geneSymbols)}, BPPARAM=param))
   }
