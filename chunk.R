library(TReNA)
library(doParallel)
library(RPostgreSQL)
source("/scratch/github/TReNA/inst/utils/createGenomeScaleModel.R")
source("/scratch/github/TReNA/inst/utils/coryGenomeScaleModel.R")

driver <- PostgreSQL()
host <- "localhost"
dbname <- "hg38"
user <- "trena"
password <- "trena"
genome.db <- dbConnect(driver, host=host, dbname=dbname, user=user, password=password)

gtex.pc.fib2 <- readRDS("/scratch/data/gtex.pc.fib2")
#fp.10 <- readRDS("/scratch/data/fp.10.rds")

list.1 <- head(rownames(gtex.pc.fib2), 1000)
system.time(fp.10 <- stinkyFeet(gtex.pc.fib2, list.1, 
   "postgres://localhost/hg38", 
   "postgres://localhost/skin_hint", 
   5000, 
   5000, 
   num.cores = 64,
   extraArgs = list("solver.list"=c("lasso", "ridge", "randomForest", "sqrtlasso", "lassopv", "pearson", "spearman"),
   "sqrtlasso"=list(num.cores=1))))

system.time(test.fib.trn <- createSpecialModel(gtex.pc.fib2, fp.10, 
 num.cores = 20,
 extraArgs = list("solver.list"=c("lasso", "ridge", "randomForest", "sqrtlasso", "lassopv", "pearson", "spearman"),
 "sqrtlasso"=list(num.cores=5))))

saveRDS(test.fib.trn, "test.fib.trn.rds")
