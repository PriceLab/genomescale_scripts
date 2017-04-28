library(TReNA)
library(doParallel)
library(RPostgreSQL)
source("/scratch/github/TReNA/inst/utils/createGenomeScaleModel.R")
source("/scratch/github/TReNA/inst/utils/coryGenomeScaleModel.R")

args <- commandArgs()
start <- as.integer(args[4])
stop <- as.integer(args[5])
name <- args[6]

driver <- PostgreSQL()
host <- "localhost"
dbname <- "hg38"
user <- "trena"
password <- "trena"
genome.db <- dbConnect(driver, host=host, dbname=dbname, user=user, password=password)

mtx <- readRDS("/scratch/data/mtx/mayo.tcx.rds")

list.1 <- rownames(mtx[start:stop,])
system.time(fp <- stinkyFeet(mtx, list.1,
   "postgres://localhost/hg38",
   "postgres://localhost/brain_hint",
   5000,
   5000,
   num.cores = 60,
   extraArgs = list("solver.list"=c("lasso", "ridge", "randomForest", "sqrtlasso", "lassopv", "pearson", "spearman"),
   "sqrtlasso"=list(num.cores=1))))

system.time(trn <- createSpecialModel(mtx, fp,
 num.cores = 20,
 extraArgs = list("solver.list"=c("lasso", "ridge", "randomForest", "sqrtlasso", "lassopv", "pearson", "spearman"),
 "sqrtlasso"=list(num.cores=5))))

saveRDS(trn, name)
