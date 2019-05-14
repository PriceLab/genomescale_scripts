# tryCatch loading files

filenames <- list.files(path = "/tmp/MODELS.cory.liver", pattern = "RData", full.names = TRUE)

# this file doesn't load
broken.file <- filenames[5463]
print(load(broken.file))

# grab a short list that includes one file that doesn't load (FAM241B)
f.names <- filenames[5460:5464]

tbls.all <- list()
for(full.path in f.names){
  if(file.info(full.path)$size==0) next
  possibleError <- tryCatch(
    {
      print(load(full.path))
    },
      # ... but if an error occurs, tell me what happened: 
      error=function(error_message) {
      message("model may be broken:", full.path)
      message("And below is the error message from R:")
      message(error_message)
      }
  )
  if(inherits(possibleError, "error")) next 
  gene.gene <- strsplit(full.path, "/", fixed=TRUE)[[1]][4]
  targetGene <- strsplit(gene.gene, ".", fixed=TRUE)[[1]][1]
  
  tbl <- results$model
  if(nrow(tbl) == 0) next
  tbl <- tbl[order(abs(tbl$pearsonCoeff), decreasing=TRUE),]
  
  tbl$targetGene <- targetGene
  tbl$rank <- seq_len(nrow(tbl))
  printf("model for %s: %d rows.orig, %d rows.trimmed, %d cols", targetGene, nrow(results$model), nrow(tbl), ncol(tbl))
  tbls.all[[targetGene]] <- tbl
}

