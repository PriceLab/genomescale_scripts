library(BiocParallel)
#----------------------------------------------------------------------------------------------------
goEnrichment <- function(geneSymbols)
{
  require("GO.db")
  require("GOstats")

  source("https://raw.githubusercontent.com/PriceLab/genomescale_scripts/master/symToGeneID.R");
  symbol.entrez.map <- assignGeneIDs(geneSymbols)
  geneIDs <- unlist(symbol.entrez.map$mapped, use.names=FALSE)


  gene.universe = character(0)

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
demo <- function()
{
   geneSymbols <- get(load("test.geneSymbols.RData"))
   length(geneSymbols)

   maxWorkers <- 40
   param=MulticoreParam(workers=maxWorkers, stop.on.error=FALSE)

   f <- function(i){   # i will be ignored; bplappy wants it
      random.indices <- sort(sample(1:length(geneSymbols), 10))
      random.geneSymbols <- geneSymbols[random.indices]
      return(goEnrichment(random.geneSymbols))
      }

   reps <- 400
   print(system.time(results <- bplapply(1:reps, f, BPPARAM=param)))

   invisible(results)

} # demo
#----------------------------------------------------------------------------------------------------

