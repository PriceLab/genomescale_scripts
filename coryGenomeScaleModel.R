library(trena)
library(BiocParallel)
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
stinkyFeet <- function(mtx.assay, gene.list, genome.db.uri, project.db.uri,
                                   size.upstream=1000, size.downstream=1000, num.cores = NULL,
                                   extraArgs = list()){

    footprint.filter <- FootprintFilter(mtx.assay = mtx.assay)

    # Setup the parallel structure with a default of half the cores
    if(is.null(num.cores)){
        num.cores <- detectCores()/2}
    cl <- makeForkCluster(nnodes = num.cores)
    registerDoParallel(cl)

    full.result.list <- foreach(i = 1:length(gene.list), .packages='TReNA', .errorhandling="pass") %dopar% {
	
	Sys.sleep(runif(1, 0, 10))
        # Designate the target gene and grab the tfs
        target.gene <- gene.list[[i]]
        out.list <- try(getCandidates(footprint.filter,extraArgs = list(
                                                  "target.gene" = target.gene,
                                                  "genome.db.uri" = genome.db.uri,
                                                  "project.db.uri" = project.db.uri,
                                                  "size.upstream" = size.upstream,
                                                  "size.downstream" = size.downstream)),
                        silent = TRUE)
	return(out.list$tfs)
	}
        # Solve the trena problem using the supplied values and the ensemble solver

    # Stop the cluster
    stopCluster(cl)

    # Name the list after the genes supplied
    names(full.result.list) <- gene.list
    return(full.result.list)

} # stinkyFeet
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
