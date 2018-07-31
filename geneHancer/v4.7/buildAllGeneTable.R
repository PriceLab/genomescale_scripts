# get all the genes on chromsome 2
# build a multi-gene table with all of the genehancer data
source("getAllEnhancers.R")
library(org.Hs.eg.db)
library(BiocParallel)
tbl.chrom <- as.data.frame(org.Hs.egCHR)
genes.chrom2 <- subset(tbl.chrom, chromosome=="2")$gene_id   # 4245
genes.all <- unique(tbl.chrom$gene_id)
length(genes.all) # 59973
keytypes(org.Hs.eg.db)
#  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"
#  [7] "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"
# [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"
# [19] "PFAM"         "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"
# [25] "UNIGENE"      "UNIPROT"
geneSyms.chrom2 <- select(org.Hs.eg.db, keys=genes.chrom2, keytype="ENTREZID", columns="SYMBOL")$SYMBOL  # 4245

length(geneSyms.chrom2)
x <- bplapply(geneSyms.chrom2, getAllEnhancers)

geneSyms.all <- unique(select(org.Hs.eg.db, keys=genes.all, keytype="ENTREZID", columns="SYMBOL")$SYMBOL)
length(geneSyms.all)  # 59905
x <- bplapply(geneSyms.all, getAllEnhancers)
length(x) # 59905
tbl.enhancers <- do.call(rbind, x)
length(unique(tbl.enhancers$geneSymbol)) # [1] 45920
save(tbl.enhancers, file="tbl.enhancers.allGenes.RData")

tbl.enhancers.chr2 <- do.call(rbind, x)
dim(tbl.enhancers.chr2)
length(unique(tbl.enhancers.chr2$geneSymbol))
save(tbl.enhancers.chr2, file="tbl.enhancer.chr2.RData")




