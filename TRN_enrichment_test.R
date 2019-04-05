# Network enrichment test to evaluate a TRN

library(GSEABase)
library(GOstats)
library(GO.db)
library(Category)
library(org.Hs.eg.db)
library(KEGG.db)
library(RUnit)

#load the whole_blood trn (GTEx)
print(load("~/Downloads/lympho-footprints-2019feb22-tbl.models.all.RData"))
whole.b <- tbl.models

#----------------------------------------------------------------------------
# pull out the top targets for a given TF (without rank)
getTF <- function(trn, geneA)
{
  temp <- subset(trn, tf.symbol == geneA)
  temp[order(temp$rank, decreasing=FALSE),]
}
#----------------------------------------------------------------------------
source("https://raw.githubusercontent.com/PriceLab/genomescale_scripts/master/symToGeneID.R"); test_assignGeneIDs()
#---------------------------------------------------------------------------------
# requires a list of genes
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
  tbl.go$genes <- unlist(geneSymbols, use.names=FALSE)
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
#---------------------------------------------------------------------------------

go.test <- function(trn, tf){
  
  target.list <- getTF(trn, tf)
  top.targets <- target.list %>% filter(rank <= 3) 
  top.targets2 <- top.targets$target.symbol
  test.go <- goEnrich(top.targets2)
  five.mean <- mean(test.go$Pvalue[1:5])
  return(five.mean)
}

go.random <- function(trn, list.length){
  
  gene.list <- unique(trn$target.symbol) %>% sample(list.length)
  test.go <- goEnrich(gene.list)
  five.mean <- mean(test.go$Pvalue[1:5])
  return(five.mean)
}
#---------------------------------------------------------------------------------
test_go.test <- function(){
  
  print(load("~/Downloads/spi1.targets.RData"))
  test.out <- go.test(spi1.targets, "SPI1")
  return(test.out)
} #1.089959e-20
