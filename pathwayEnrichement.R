library(GSEABase)
library(GOstats)
library(GO.db)
library(Category)
library(org.Hs.eg.db)
library(KEGG.db)
#---------------------------------------------------------------------------------
# kegg and GO enrichment needs entrez geneIDs.  we usually start with gene sybmols
# conversion is offered in this file:
# note that the assignGeneIDs function provide here returns a list with three
# name fields:
#     mapped: successful, unique entrezIDs found for symbols
#     multiples: quasi-successful, multiple ids found for each of these symbols
#     failures: no entrezIDs found for these symbols
source("symToGeneID.R"); test_assignGeneIDs()
#---------------------------------------------------------------------------------
# requires a list of genes
#---------------------------------------------------------------------------------
goEnrich <- function(geneSymbols)
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
    tbl.go$genes <- geneSymbols
    tbl.go

} # goEnrich
#---------------------------------------------------------------------------------
keggEnrich <- function(geneIDs)
{
   symbol.entrez.map <- assignGeneIDs(geneIDs)

   gene.universe = character(0)
   geneIDs <- unlist(symbol.entrez.map$mapped, use.names=FALSE)

   kegg.params <- new("KEGGHyperGParams", geneIds = unique(geneIDs),
                      universeGeneIds = character(0), annotation = "org.Hs.eg.db",
                      pvalueCutoff = 0.1, testDirection = "over")

   kegg.hgr  <- hyperGTest(kegg.params)

   return(kegg.hgr)

} # keggEnrich
#--------------------------------------------------------------------------------
test_goEnrich <- function()
{
   igap.ad.genes <- c("CR1", "BIN1", "CD2AP", "EPHA1", "CLU", "MS4A6A", "PICALM",
                      "ABCA7", "CD33", "HLA-DRB5", "HLA-DRB1", "PTK2B", "SORL1",
                      "SLC24A4", "RIN3", "DSG2", "INPP5D", "MEF2C", "NME8", "ZCWPW1",
                      "CELF1", "FERMT2", "CASS4", "APOE", "TOMM40")

   #x <- assignGeneIDs(igap.ad.genes)
   #lapply(x, length)   # mapped: 25   failures: 0   muplitples: 0
   #geneIDs <- unlist(x$mapped, use.names=FALSE)

   tbl.go <- goEnrich(head(igap.ad.genes))


} # test_goEnrich
#--------------------------------------------------------------------------------

