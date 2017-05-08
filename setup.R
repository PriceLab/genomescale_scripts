source("https://bioconductor.org/biocLite.R")
library(BiocInstaller)
biocLite(c("glmnet", "GenomicRanges", "RSQLite", "lassopv", "randomForest", "flare", "vbsr", "foreach", "doParallel", "RPostgreSQL", "RMySQL", "DBI", "RUnit", "BiocParallel", "org.Hs.eg.db", "Biostrings", "rtracklayer", "dplyr"), ask=FALSE)
