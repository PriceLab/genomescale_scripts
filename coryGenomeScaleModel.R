library(trena)
library(BiocParallel)
library(RPostgreSQL)
library(dplyr)
#----------------------------------------------------------------------------------------------------
# Bring in the TF-motif mapping
motifsgenes <- readRDS("./2017_10_26_Motif_TF_Map.RDS")
#----------------------------------------------------------------------------------------------------
createGenomeScaleModel <- function(mtx.assay,
                                   gene.list,
                                   genome.db.uri,
                                   project.db.uri,
                                   size.upstream=1000,
                                   size.downstream=1000,
                                   num.cores = NULL,
                                   nCores.sqrt = 4,
                                   solverNames){
    
    lapply(dbListConnections(dbDriver(drv="PostgreSQL")), dbDisconnect)

    # Setup the parallel structure with a default of half the cores
    if(is.null(num.cores)){
        num.cores <- detectCores()/2}

    # Use BiocParallel
    register(MulticoreParam(workers = num.cores,
                            stop.on.error = FALSE,
                            log = TRUE),
             default = TRUE)

    # Make the model-creation into a function
    createGeneModel <- function(target.gene, mtx.assay, genome.db.uri, project.db.uri,
                                size.upstream, size.downstream, solverNames){

        # Create the footprint filter and get candidates with it
        footprint.filter <- FootprintFilter(genomeDB = genome.db.uri,
                                            footprintDB = project.db.uri,                                          
                                            geneCenteredSpec = list(targetGene = target.gene,
                                                                    tssUpstream = size.upstream,
                                                                    tssDownstream = size.downstream),
                                            regionsSpec = list())        
        out.list <- try(getCandidates(footprint.filter),silent = TRUE)

        # Solve the trena problem using the supplied values and the ensemble solver
        if(!(class(out.list) == "try-error")){
            if(length(out.list$tfs) > 0){
                trena <- EnsembleSolver(mtx.assay, 
					targetGene = target.gene,
					candidateRegulators = out.list$tfs, 
                                        solverNames = solverNames,                                        
                                        nCores.sqrt = nCores.sqrt)                
                return(solve(trena))
            }            
            else{return(NULL)}            
        }
        else{return(NULL)}                
    }

    # Run the function for the gene list using bplapply
    result <- bplapply(gene.list, createGeneModel,
                       mtx.assay = mtx.assay,
                       genome.db.uri = genome.db.uri,
                       project.db.uri = project.db.uri,                       
                       size.upstream = size.upstream,                       
                       size.downstream = size.downstream,                       
                       solverNames = solverNames)
    return(result)

} # createGenomeScaleModel
#----------------------------------------------------------------------------------------------------
# Note: Run this on a dataframe of regions, including gene names (geneSymbol column)
getTfsFromDb <- function(regions, genome.db.uri, project.db.uri,
                       size.upstream=5000, size.downstream=5000, num.cores = 8){

    # Setup the parallel structure with a default of half the cores
#    if(is.null(num.cores)){
#        num.cores <- detectCores()/2}

    # Use BiocParallel    
    register(MulticoreParam(workers = num.cores,
#    register(SerialParam(    

                            stop.on.error = FALSE,                            
                            log = TRUE),
             default = TRUE)

    # Transform the given regions into a list of region dataframes
    regions.wo.genes <- select(regions, -geneSymbol)

    # Make the dataframe into a list of dataframes
    dfToList <- function(regions){
        df.list <- list()
        for(i in 1:floor(nrow(regions)/10)){
            idx1 <- 10*i-9
            idx2 <- 10*i
            df.list[[i]] <- regions[idx1:idx2,]
        }
        
        if(nrow(regions) %% 10 != 0){
            i <- floor(nrow(regions)/10)
            idx1 <- 10*i+1
            idx2 <- nrow(regions)
            df.list[[i+1]] <- regions[idx1:idx2,]
        }
        
        return(df.list)
    }    

    regions.list <- dfToList(regions)
    
    # Function to convert motifs to tfs
    convertMotifsToTfs <- function(motifs){

        # Catch footprints that don't exist
        if(is.character(motifs)) return(NA)
        tf.df <- motifsgenes %>%
            filter(motif %in% motifs$motifName)
        return(unique(tf.df$tf))                
    }

    selectOrNA <- function(output){
        # If it's a dataframe, return the motifName column
        if(is.character(output)){
            return(output)
        } else if(nrow(output) == 0){
            return("No footprints found")}        
        return(select(output, motifName))
    }
                
    findGeneFootprints <- function(regions, genome.db.uri, project.db.uri){
       
        # Create the footprint filter from the target gene
        footprint.filter <- try(FootprintFilter(genomeDB = genome.db.uri,
                                                footprintDB = project.db.uri,
                                                regions = regions),                                
                                silent = TRUE)        

        # Only grab candidates if the filter is valid
        if(class(footprint.filter) == "FootprintFilter"){
            out.list <- getCandidates(footprint.filter)

            # Catch empty lists
            if(length(out.list) == 0) return(character(0))

            # Only return TFs if candidate grab is not null
            if(class(out.list) != "NULL"){
                # Use a semi join to grab the correct tfs
                motif.list <- lapply(out.list, selectOrNA)
                tf.list <- lapply(motif.list, convertMotifsToTfs)

                return(tf.list)

            } else {
                return("No Candidates Found")
                }                
        } else{            
            return(footprint.filter[1])
            }
    }    

    full.result.list <- bplapply(regions.list, findGeneFootprints,
                                 genome.db.uri = genome.db.uri,
                                 project.db.uri = project.db.uri)
    
    # Un-nest and Name the list after the genes supplied
    full.result.list <- unlist(full.result.list, recursive = FALSE)    
    names(full.result.list) <- regions$geneSymbol

    # Remove any where the content is wrong
    no.fp <- which(!(sapply(full.result.list, is.character)))
    full.result.list[no.fp] <- NULL
    
    return(full.result.list)

} # getTfsFromDb
#------------------------------------------------------------------------------------------------------
createSpecialModel <- function(mtx.assay, gene.list, num.cores = NULL,
                                   extraArgs = list()){

    trena <- TReNA(mtx.assay, solver = "ensemble")

    #lapply(dbListConnections(dbDriver(drv="PostgreSQL")), dbDisconnect)

    # Setup the parallel structure with a default of half the cores
    if(is.null(num.cores)){
        num.cores <- detectCores()/2}

    
    cl <- makePSOCKcluster(num.cores)
    registerDoParallel(cl)

    full.result.list <- foreach(i = 1:length(names(gene.list)), .packages='TReNA', .errorhandling="pass") %dopar% {

        # Designate the target gene and grab the tfs
        target.gene <- names(gene.list)[[i]]

        # Solve the trena problem using the supplied values and the ensemble solver

        if(!(class(gene.list[[target.gene]]) == "try-error")){
            if(length(gene.list[[target.gene]]$tfs) > 0){

                solve(trena, target.gene, gene.list[[target.gene]]$tfs, extraArgs = extraArgs)}

            else{NULL}


        }
        else{NULL}
}
    # Stop the cluster
    stopCluster(cl)

    # Name the list after the genes supplied
    names(full.result.list) <- names(gene.list)
    return(full.result.list)

} # createSpecialModel
#----------------------------------------------------------------------------------------------------
getTfsFromAllDbs <- function(mtx.assay, gene.list, genome.db.uri, project.list,
                             size.upstream=1000, size.downstream=1000, num.cores = NULL)
{
    footprint.filter <- FootprintFilter(mtx.assay = mtx.assay)

    # Setup the parallel structure with a default of half the cores
    if(is.null(num.cores)){
        num.cores <- detectCores() - 1}
    cl <- makeForkCluster(num.cores)
    registerDoParallel(cl)


    # Pass the appropriate variables
#    clusterExport(cl, varlist = c("footprint.filter","gene.list",
#                                  "genome.db.uri","project.list","size.upstream",
#                                  "size.downstream"))
    result.list <- foreach(i = 1:length(gene.list)) %dopar% {
      #  1}

	#Sys.sleep(runif(1, 0, 10))
        # Designate the target gene and grab the tfs only from each of the 4 databases
        my.target <- gene.list[[i]]
        all.tfs <- character()

        # Loop through the list of project dbs and grab tfs from each
        for(project in project.list){
            out.list <- try(getCandidates(footprint.filter,extraArgs = list(
                                                               "target.gene" = my.target,
                                                               "genome.db.uri" = genome.db.uri,                                                              
                                                               "project.db.uri" = project,
                                                               "size.upstream" = size.upstream,
                                                               "size.downstream" = size.downstream)),                                                                      
                            silent = TRUE)
            # Add to the list only if it has tfs
            if(!(class(out.list) == "try-error")){
                if(length(out.list$tfs) > 0){
                    all.tfs <- c(all.tfs,out.list$tfs)
                }
            }            
        }
        # Return the union
                return(unique(all.tfs))
    }
    
    # Stop the cluster
    stopCluster(cl)

    # Name the list after the genes supplied
    names(result.list) <- gene.list
    return(result.list)
} # getTfsFromAllDbs
#----------------------------------------------------------------------------------------------------
createAverageModel <- function(mtx.assay, gene.list, num.cores = NULL,
                                   extraArgs = list()){

    trena <- TReNA(mtx.assay, solver = "ensemble")

    #lapply(dbListConnections(dbDriver(drv="PostgreSQL")), dbDisconnect)

    # Setup the parallel structure with a default of half the cores
    if(is.null(num.cores)){
        num.cores <- detectCores() - 1}
    cl <- makePSOCKcluster(num.cores)
    registerDoParallel(cl)

    full.result.list <- foreach(i = 1:length(names(gene.list)), .packages='TReNA', .errorhandling="pass") %dopar% {

        # Designate the target gene and grab the tfs
        target.gene <- names(gene.list)[[i]]

        # Solve the trena problem using the supplied values and the ensemble solver
        if(!(class(gene.list[[target.gene]]) == "try-error")){
            if(length(gene.list[[target.gene]]) > 0){

                solve(trena, target.gene, gene.list[[target.gene]], extraArgs = extraArgs)}

            else{NULL}
        }
        else{NULL}
}
    # Stop the cluster
    stopCluster(cl)

    # Name the list after the genes supplied
    names(full.result.list) <- names(gene.list)
    return(full.result.list)

} # createAverageModel
#----------------------------------------------------------------------------------------------------
createModelFromGeneList <- function(mtx.assay, gene.list, num.cores = NULL,
                                    solverList = c("lasso","ridge"),
                                    nCores.sqrt = 2){

    # Remove genes from the list that don't have any TFs
    rm.idx <- which(sapply(gene.list,length) == 1)
    gene.list[rm.idx] <- NULL
    
    # Create parallel structure w/ BiocParallel
    register(MulticoreParam(workers = num.cores,                            
                            stop.on.error = FALSE,                            
                            log = TRUE),             
             default = TRUE)


    # Create a function that:
    # 1) Takes a Named List (name = target.gene, list = regulators)
    # 2) Creates an ensemble solver with the prescribed solvers
    # 3) Solves the solver

    buildAndSolveForGene <- function(idx,gene.list, mtx.assay, solverList, nCores.sqrt){

        # Build the ensemble solver
        e.solver <- EnsembleSolver(mtx.assay = mtx.assay,
                                   targetGene = names(gene.list)[idx],
                                   candidateRegulators = gene.list[[idx]],
                                   solverNames = solverList,
                                   nCores.sqrt = nCores.sqrt)
        # Solve the ensemble solver
        return(run(e.solver))
        }

    full.result.list <- bptry(bplapply(1:length(gene.list), buildAndSolveForGene,
                                       gene.list = gene.list,
                                       mtx.assay = mtx.assay,
                                       solverList = solverList,
                                       nCores.sqrt = nCores.sqrt
                                       )
                              )
    # Name the list after the genes supplied
    names(full.result.list) <- names(gene.list)
    return(full.result.list)

} # createModelFromGeneList
#----------------------------------------------------------------------------------------------------
getTfsFromSampleIDs <- function(gene.list, sampleIDs, genome.db.uri, project.db.uri,
                                size.upstream=1000, size.downstream=1000, num.cores = 8){

    # Setup the parallel structure with a default of half the cores
    if(is.null(num.cores)){
        num.cores <- detectCores()/2}

    # Use BiocParallel    
    register(MulticoreParam(workers = num.cores,
    #register(SerialParam(    

                            stop.on.error = FALSE,                            
                            log = TRUE),
             default = TRUE)

    findGeneFootprints <- function(target.gene, genome.db.uri, project.db.uri,
                                   size.upstream, size.downstream, sampleIDs){

        # Create the footprint filter from the target gene
        footprint.filter <- try(FootprintFilter(genomeDB = genome.db.uri,
                                                footprintDB = project.db.uri,
                                                geneCenteredSpec = list(targetGene = target.gene,
                                                                        tssUpstream = size.upstream,
                                                                        tssDownstream = size.downstream),
                                                regionsSpec = list()),                                
                                silent = TRUE)

        # Only grab candidates if the filter is valid
        if(class(footprint.filter) == "FootprintFilter"){
            out.list <- getCandidates(footprint.filter)

            # Only return TFs if candidate grab is not null
            if(class(out.list) != "NULL"){

                # Filter out only the desired sampleIDs
                out.list$tbl <- filter(out.list$tbl, sample_id %in% sampleIDs)
                out.list$tfs <- unique(out.list$tbl$tf)                
                
                return(out.list$tfs)
            } else {
                return("No Candidates Found")
                }                
        } else{            
            return(footprint.filter[1])
            }
    }    

    full.result.list <- bplapply(gene.list, findGeneFootprints,
                                 genome.db.uri = genome.db.uri,
                                 project.db.uri = project.db.uri,
                                 size.upstream = size.upstream,
                                 size.downstream = size.downstream,
                                 sampleIDs = sampleIDs)
    
    # Name the list after the genes supplied
    names(full.result.list) <- gene.list
    return(full.result.list)

} # getTfsFromSampleIDs
#------------------------------------------------------------------------------------------------------
getTfsFromSampleIDsMultiDB <- function(gene.list, sampleIDs, genome.db.uri, projectList,
                                size.upstream=1000, size.downstream=1000, num.cores = 8){

    # Setup the parallel structure with a default of half the cores
    if(is.null(num.cores)){
        num.cores <- detectCores()/2}

    # Use BiocParallel    
    register(MulticoreParam(workers = num.cores,
    #register(SerialParam(    

                            stop.on.error = FALSE,                            
                            log = TRUE),
             default = TRUE)

    findGeneFootprints <- function(target.gene, genome.db.uri, project.db.uri,
                                   size.upstream, size.downstream, sampleIDs){
        
        # Create the footprint filter from the target gene
        footprint.filter <- try(FootprintFilter(genomeDB = genome.db.uri,                                                
                                                footprintDB = project.db.uri,                                             
                                                geneCenteredSpec = list(targetGene = target.gene,
                                                                        tssUpstream = size.upstream,
                                                                        tssDownstream = size.downstream),
                                                regionsSpec = list()),                                
                                silent = TRUE)                    
        # Only grab candidates if the filter is valid        
        if(class(footprint.filter) == "FootprintFilter"){            
            out.list <- getCandidates(footprint.filter)                            
            # Only return TFs if candidate grab is not null            
            if(class(out.list) != "NULL"){                                    
                # Filter out only the desired sampleIDs                
                out.list$tbl <- filter(out.list$tbl, sample_id %in% sampleIDs)                
                out.list$tfs <- unique(out.list$tbl$tf)                                                    
                return(out.list$tfs)                
            } else {                
                return("No Candidates Found")                
            }            
        } else{            
            return("Cannot create filter")            
        }        
    }

    # Define a function that loops through a list and accumulates TF lists
    combineTFsFromDBs <- function(target.gene, genome.db.uri, projectList,
                                  size.upstream, size.downstream, sampleIDs){

        # Create an empty vector
        all.tfs <- character(0)

        # Find Footprints from each DB and add to list
        for(project.db.uri in projectList){
            new.tfs <- findGeneFootprints(target.gene, genome.db.uri, project.db.uri,
                                          size.upstream, size.downstream, sampleIDs)
            all.tfs <- union(all.tfs, new.tfs)
        }
        # Return the full list
        return(all.tfs)
    }        
    
    # This part should remain the same
    full.result.list <- bplapply(gene.list, combineTFsFromDBs,
                                 genome.db.uri = genome.db.uri,
                                 projectList = projectList,
                                 size.upstream = size.upstream,
                                 size.downstream = size.downstream,
                                 sampleIDs = sampleIDs)
    
    # Name the list after the genes supplied
    names(full.result.list) <- gene.list
    return(full.result.list)

} # getTfsFromSampleIDsMultiDB
#----------------------------------------------------------------------------------------------------
# Note: Run this on a dataframe of regions, including gene names (geneSymbol column)
getTfsFromMultiDB <- function(regions, genome.db.uri, projectList,num.cores = 8){
        
    # Make the dataframe into a list of dataframes
    dfToList <- function(regions){
        df.list <- list()
        for(i in 1:floor(nrow(regions)/10)){
            idx1 <- 10*i-9
            idx2 <- 10*i
            df.list[[i]] <- regions[idx1:idx2,]
        }
        
        if(nrow(regions) %% 10 != 0){
            i <- floor(nrow(regions)/10)
            idx1 <- 10*i+1
            idx2 <- nrow(regions)
            df.list[[i+1]] <- regions[idx1:idx2,]
        }
        
        return(df.list)
    }    

    regions.list <- dfToList(regions)
    
    # Function to convert motifs to tfs
    convertMotifsToTfs <- function(motifs){

        # Catch footprints that don't exist
        if(is.character(motifs)) return(NA)
        tf.df <- motifsgenes %>%
            filter(motif %in% motifs$motifName)
        return(unique(tf.df$tf))                
    }

    selectOrNA <- function(output){
        # If it's a dataframe, return the motifName column
        if(is.character(output)){
            return(output)
        } else if(nrow(output) == 0){
            return("No footprints found")}        
        return(select(output, motifName))
    }
                
    findGeneFootprints <- function(regions, genome.db.uri, project.db.uri){
       
        # Create the footprint filter from the target gene
        footprint.filter <- try(FootprintFilter(genomeDB = genome.db.uri,
                                                footprintDB = project.db.uri,
                                                regions = regions),                                
                                silent = TRUE)        

        # Only grab candidates if the filter is valid
        if(class(footprint.filter) == "FootprintFilter"){
            out.list <- getCandidates(footprint.filter)

            # Catch empty lists
            if(length(out.list) == 0) return(character(0))

            # Only return TFs if candidate grab is not null
            if(class(out.list) != "NULL"){
                # Use a semi join to grab the correct tfs
                motif.list <- lapply(out.list, selectOrNA)
                tf.list <- lapply(motif.list, convertMotifsToTfs)

                return(tf.list)

            } else {
                return("No Candidates Found")
                }                
        } else{            
            return(footprint.filter[1])
            }
    }    
    
    # Define a function that loops through a list and accumulates TF lists
    combineTFsFromDBs <- function(regions, genome.db.uri, projectList){

        # Take in the regions DF with gene symbol and pull it off
        regions.wo.genes <- select(regions, -geneSymbol)        

        # Find the first set of footprints
        all.tfs <- findGeneFootprints(regions.wo.genes,
                                      genome.db.uri,
                                      projectList[1])

        # Name the list after the genes supplied        
        names(all.tfs) <- regions$geneSymbol

        # Find Footprints from each DB and add to list
        for(i in 1:length( projectList)){
            new.tfs <- findGeneFootprints(regions.wo.genes,
                                          genome.db.uri,
                                          projectList[i])

            # Name and un-nest the TFs as before          
            names(new.tfs) <- regions$geneSymbol

            # Consolidate the 2 lists
            keys <- names(all.tfs)
            all.tfs <- setNames(mapply(union,
                                       all.tfs[keys],
                                       new.tfs[keys]),
                                keys)
        }
        # Return the full list
        return(all.tfs)
    }        

    # Use BiocParallel    
    register(MulticoreParam(workers = num.cores,
                            stop.on.error = FALSE,                            
                            log = TRUE),
             default = TRUE)

    full.result.list <- bplapply(regions.list, combineTFsFromDBs,
                                 genome.db.uri = genome.db.uri,
                                 projectList = projectList)

    # Un-nest the list
    full.result.list <- unlist(full.result.list, recursive = FALSE)
    
    # Remove any where the content is wrong
    no.fp <- which(!(sapply(full.result.list, is.character)))
    full.result.list[no.fp] <- NULL
    return(full.result.list)
    
} # getTfsFromSampleIDsMultiDB
#----------------------------------------------------------------------------------------------------
# Example Cory Script
# Assume my.mtx is the matrix, hg38 is the genome.db, brain is the tissue, shoulder is 5000

# Also assume this has 128 cores!!
# Step 1: Get all the genes

testRun <- function(my.mtx){

    all.genes <- getTfsFromMultiDB(rownames(my.mtx),
                                   genome.db.uri = "postgres://localhost/hg38",
                                   projectList = c("postgres://localhost/brain_hint_20",
                                                   "postgres://localhost/brain_hint_16",
                                                   "postgres://localhost/brain_wellington_20",
                                                   "postgres://localhost/brain_wellington_16"),
                                   size.upstream = 5000,
                                   size.downstream = 5000,
                                   num.cores = 100)
    
    # Step 2: Use all the genes to make ALL the models
    all.models <- createModelFromGeneList(my.mtx, all.genes, num.cores = 30,
                                          solverList = c("lasso","ridge","pearson",
                                                         "spearman","randomforest",
                                                         "lassopv","sqrtlasso"),
                                          nCores.sqrt = 4)
}
#----------------------------------------------------------------------------------------------------
getProxProbesPromoter <- function(probeIDs,
                   tssUpstream = 5000,
                   tssDownstream = 5000){

              # Switch the name of the database and filter we use
              db.name <-  "hsapiens_gene_ensembl"
              
              filter.name <- "illumina_humanht_12_v4"

              my.mart <- biomaRt::useMart(biomart="ensembl", dataset= db.name)

              tbl.geneInfo <- biomaRt::getBM(attributes=c("chromosome_name",
                                                          "transcription_start_site",
                                                          "transcript_tsl",
                                                          "hgnc_symbol",
                                                          filter.name),
                                             filters=filter.name, value=probeIDs, mart=my.mart)

              if(nrow(tbl.geneInfo) == 0)
                  return(NA)

              # Sort by hgnc_symbol and transcript_tsl, then pull the first entry for each gene
              tbl.geneInfo <- tbl.geneInfo[order(tbl.geneInfo[[filter.name]],
                                                 tbl.geneInfo$transcript_tsl),]
              tbl.geneInfo <- tbl.geneInfo[match(unique(tbl.geneInfo[[filter.name]]),
                                                 tbl.geneInfo[[filter.name]]),]

              # remove contigs and check to make sure it's just 1 chromosome
              tbl.geneInfo <- subset(tbl.geneInfo, chromosome_name %in% c(1:22, "X", "Y", "MT"))
              chrom <- sprintf("chr%s", tbl.geneInfo$chromosome_name)

              tss <- tbl.geneInfo$transcription_start_site
              start.loc <- tss - tssDownstream
              end.loc   <- tss + tssUpstream

              temp <- data.frame(geneSymbol=tbl.geneInfo$hgnc_symbol,
                                 chrom=chrom,
                                 start=start.loc,
                                 end=end.loc,
                                 stringsAsFactors=FALSE)
                                 
              return (temp[!(duplicated(temp$geneSymbol)),])

          }

#----------------------------------------------------------------------------------------------------
# How to call it; a sample function
sampleCall <- function(regions){

    # Assume we've got a set of regions...
    genome.db.uri <- "postgres://localhost/hg38"
    projectList <- c("postgres://localhost/brain_hint_20",
                     "postgres://localhost/brain_hint_16",
                     "postgres://localhost/brain_wellington_20",
                     "postgres://localhost/brain_wellington_16")

    # Call using 30 cores
    all.candidates <- getTfsFromMultiDB(regions, genome.db.uri, projectList, 30)

} #sampleCall

# For Cory: assuming you've called your file "my.regions"
# my.stuff <- sampleCall(my.regions)

    


 
