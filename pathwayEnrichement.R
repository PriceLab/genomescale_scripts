library(GSEABase)
library(GOstats)
library(GO.db)
library(Category)
library(org.Hs.eg.db)
library(KEGG.db)

# requires a list of genes
#---------------------------------------------------------------------------------
goEnrich <- function(genes){

symbol.entrez.map <- assignGeneIDs(genes)

gene.universe = character(0)
geneIDs <- unlist(symbol.entrez.map$mapped, use.names=FALSE)

go.params <- new("GOHyperGParams", geneIds=unique(geneIDs),
                 universeGeneIds = gene.universe, annotation = "org.Hs.eg.db",
                 ontology = 'BP', pvalueCutoff = 0.05, conditional = FALSE,
                 testDirection = "over")

go.bp.hgr <- hyperGTest(go.params)

return(go.bp.hgr)
}
#---------------------------------------------------------------------------------
keggEnrich <- function(genes){

symbol.entrez.map <- assignGeneIDs(genes)

gene.universe = character(0)
geneIDs <- unlist(symbol.entrez.map$mapped, use.names=FALSE)

kegg.params <- new("KEGGHyperGParams", geneIds = unique(geneIDs),
                   universeGeneIds = character(0), annotation = "org.Hs.eg.db",
                   pvalueCutoff = 0.1, testDirection = "over")

kegg.hgr  <- hyperGTest(kegg.params)

return(kegg.hgr)
}


