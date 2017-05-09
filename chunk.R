library(TReNA)
library(doParallel)
library(RPostgreSQL)
source("/scratch/github/TReNA/inst/utils/createGenomeScaleModel.R")
source("/scratch/github/TReNA/inst/utils/coryGenomeScaleModel.R")

args <- commandArgs()
start <- as.integer(args[4])
stop <- as.integer(args[5])
name <- args[6]

mtx <- as.matrix(readRDS("/scratch/data/mayo.cer.rds"))

gene.list <- rownames(mtx[start:stop,])

genome.db.uri <- "postgres://localhost/hg38"
project.list <- c("postgres://localhost/brain_hint_16",
    		"postgres://localhost/brain_hint_20",
            	"postgres://localhost/brain_wellington_16",
                "postgres://localhost/brain_wellington_20")
size.downstream <- 5000
size.upstream <- 5000   					
num.cores <- 60
system.time(fp <- getTfsFromAllDbs(mtx, gene.list, genome.db.uri, project.list, size.upstream, size.downstream, num.cores))

fp <- fp[lapply(fp,length)>0]

system.time(trn <- createAverageModel(mtx, fp, 
 num.cores = 20,
 extraArgs = list("solver.list"=c("lasso", "ridge", "randomForest", "sqrtlasso", "lassopv", "pearson", "spearman"),
 "sqrtlasso"=list(num.cores=5))))

saveRDS(trn, name)
