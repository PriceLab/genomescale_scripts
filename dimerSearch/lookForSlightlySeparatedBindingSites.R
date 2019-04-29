library(MotifDb)
source("~/github/trenaShinyApps/utils/plotMotifs.R")
dimers <- unique(toupper(grep("::", mcols(MotifDb)$geneSymbol, fixed=TRUE, value=TRUE)))
length(dimers)
tbl.atf <- get(load("~/github/genomescale_scripts/dimerSearch/tbl.associatedTFs.combined.RData"))
subset(tbl.atf, tf.1=="TAL1" & tf.2=="TCF3")  #  TAL1 TCF3     2
pfms <- query(MotifDb, c("jaspar2018", "hsapiens"), c("TAL1", "TCF3"))
plotMotifs(pfms)

for(i in 9:length(dimers)){
   dimer <- dimers[i]
   tokens <- strsplit(dimer, "::", fixed=TRUE)[[1]]
   a <- tokens[1]
   b <- tokens[2]
   hits.ab <- nrow(subset(tbl.atf, tf.1==a & tf.2==b))
   hits.ba <- nrow(subset(tbl.atf, tf.1==b & tf.2==a))
   hits.both <- hits.ab + hits.ba
   if(hits.both > 0){
      printf("%s: %d %d", dimer, hits.ab, hits.ba)
      pfms <- query(MotifDb, c("hsapiens"), sprintf("%s::%s", a, b))
      if(length(pfms) > 0)
         plotMotifs(pfms)
      browser()
      xyz <- 99
      }
   } # for i


# jaspar_2018-stat1::stat2-Ma0517.1 might be two motifs glued together - no!  seem to have same motif
# as might be jaspar2018-ppara::rxra-ma1148.1 - no.  also seem to have same motif
