# get all the genes on chromsome 2
# build a multi-gene table with all of the genehancer data
source("getAllEnhancers.R")
library(org.Hs.eg.db)
tbl.chrom <- as.data.frame(org.Hs.egCHR)
genes.chrom2 <- subset(tbl.chrom, chromosome=="2")$gene_id   # 4245
keytypes(org.Hs.eg.db)
#  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"
#  [7] "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"
# [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"
# [19] "PFAM"         "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"
# [25] "UNIGENE"      "UNIPROT"
geneSyms.chrom2 <- select(org.Hs.eg.db, keys=genes.chrom2, keytype="ENTREZID", columns="SYMBOL")$SYMBOL  # 4245
x <- lapply(geneSyms.chrom2, getAllEnhancers)
head(genes.chrom2)



