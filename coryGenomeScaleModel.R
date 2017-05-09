createGenomeScaleModel <- function(mtx.assay, gene.list, genome.db.uri, project.db.uri,
                                   size.upstream=1000, size.downstream=1000, num.cores = NULL,
                                   extraArgs = list()){

    footprint.filter <- FootprintFilter(mtx.assay = mtx.assay)
    trena <- TReNA(mtx.assay, solver = "ensemble")

    lapply(dbListConnections(dbDriver(drv="PostgreSQL")), dbDisconnect)

    # Setup the parallel structure with a default of half the cores
    if(is.null(num.cores)){
        num.cores <- detectCores()/2}  
    cl <- makeForkCluster(nnodes = num.cores)
    registerDoParallel(cl)
    
    full.result.list <- foreach(i = 1:length(gene.list), .packages='TReNA') %dopar% {
	
	# Designate the target gene and grab the tfs
        target.gene <- gene.list[[i]]        
        out.list <- try(getCandidates(footprint.filter,extraArgs = list(
                                                  "target.gene" = target.gene,
                                                  "genome.db.uri" = genome.db.uri,
                                                  "project.db.uri" = project.db.uri,
                                                  "size.upstream" = size.upstream,                                          
                                                  "size.downstream" = size.downstream)),
                        silent = TRUE)
	
        # Solve the trena problem using the supplied values and the ensemble solver

        if(!(class(out.list) == "try-error")){
            if(length(out.list$tfs) > 0){
                
                solve(trena, target.gene, out.list$tfs, extraArgs = extraArgs)}
            
            else{NULL}

            
        }
        else{NULL}
}
    # Stop the cluster
    stopCluster(cl)
    
    # Name the list after the genes supplied
    names(full.result.list) <- gene.list
    return(full.result.list)

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
# Commands to create variables for testing above
library(TReNA);library(doParallel)
mtx.assay <- as.matrix(readRDS("/scratch/data/mayo.tcx.rds"))
gene.list <- head(rownames(mtx.assay))
genome.db.uri <- "postgres://localhost/hg38"
project.list <- c("postgres://localhost/brain_hint_16",
                  "postgres://localhost/brain_hint_20")
size.upstream <- 5000; size.downstream <- 5000; num.cores <-10

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
