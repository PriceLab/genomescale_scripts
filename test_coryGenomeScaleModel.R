library(TReNA)
library(doParallel)
library(RUnit)
source("./coryGenomeScaleModel.R")

test_getTfsFromAllDbs <- function(){

    mtx.assay <- as.matrix(readRDS("/scratch/data/mayo.tcx.rds"))    
    gene.list <- head(rownames(mtx.assay),10)
    genome.db.uri <- "postgres://localhost/hg38"    
    project.list <- c("postgres://localhost/brain_hint_16",                      
                      "postgres://localhost/brain_hint_20")    
    size.upstream <- 5000; size.downstream <- 5000; num.cores <-10

    # Use the getTfsFromAllDbs function
    all.results <- getTfsFromAllDbs(mtx.assay, gene.list, genome.db.uri, project.list,
                                    size.upstream, size.downstream, num.cores)
    
    # Check that lengths are correct    
    checkEquals(length(all.results[[1]]),0)
    checkEquals(length(all.results[[2]]),0)
    checkEquals(length(all.results[[3]]),0)
    checkEquals(length(all.results[[4]]),360)
    checkEquals(length(all.results[[5]]),360)
    checkEquals(length(all.results[[6]]),438)

    # Check that names are correct
    checkEquals(names(all.results),gene.list)
                
} # test_getTfsFromAllDbs
