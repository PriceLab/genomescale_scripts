library(org.Hs.eg.db)
library(Category)
library(GOstats)

geneIDs <- c("4487",  "54345", "2737",  "64321", "83595", "30812", "22887", "79192", "4092", "3664")
gene.universe = character(0)

go.params <- new("GOHyperGParams", geneIds=unique(geneIDs),
                 universeGeneIds = gene.universe, annotation = "org.Hs.eg.db",
                 ontology = 'BP', pvalueCutoff = 0.05, conditional = FALSE,
                 testDirection = "over")

print(go.params@datPkg@db)
slotNames(go.params@datPkg)
go.params@datPkg@name
go.params@datPkg@db
go.params@datPkg@installed

RSQLite::dbDisconnect(dbconn(go.params@datPkg@db))

rm(go.params)

go.params2 <- new("GOHyperGParams", geneIds=unique(geneIDs),
                 universeGeneIds = gene.universe, annotation = "org.Hs.eg.db",
                 ontology = 'BP', pvalueCutoff = 0.05, conditional = FALSE,
                 testDirection = "over")

go.params2@datPkg@name
go.params2@datPkg@installed
go.params2@datPkg@db

