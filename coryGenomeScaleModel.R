library(trena)
library(BiocParallel)
library(RPostgreSQL)
library(dplyr)
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
getTfsFromDb <- function(gene.list, genome.db.uri, project.db.uri,
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
                                   size.upstream, size.downstream){

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
                                 size.downstream = size.downstream)
    
    # Name the list after the genes supplied
    names(full.result.list) <- gene.list
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
getTfsFromMultiDB <- function(gene.list, genome.db.uri, projectList,
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
                                   size.upstream, size.downstream){
        
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
                                  size.upstream, size.downstream){

        # Create an empty vector
        all.tfs <- character(0)

        # Find Footprints from each DB and add to list
        for(project.db.uri in projectList){
            new.tfs <- findGeneFootprints(target.gene, genome.db.uri, project.db.uri,
                                          size.upstream, size.downstream)
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
# Example Cory Script
# Assume my.mtx is the matrix, hg38 is the genome.db, brain is the tissue, shoulder is 5000

# Also assume this has 128 cores!!
# Step 1: Get all the genes
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
