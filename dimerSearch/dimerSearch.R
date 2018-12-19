library(RUnit)
library(gplots)   # for heatmap.2
library(PSICQUIC)
library(trenaSGM)  # if you want to use allKnownTFs()
#------------------------------------------------------------------------------------------------------------------------
if(!exists("psicquic")){
   psicquic <- PSICQUIC(test=FALSE)
   providers <- providers(psicquic)
   psicquic.id.mapper <- IDMapper("9606")
   }

#------------------------------------------------------------------------------------------------------------------------
tbl.models <- get(load("tbl.models.all.RData"))
#------------------------------------------------------------------------------------------------------------------------
# extract all models which have <tf> at or below rank.max
# determine which tfs accompany <tf> at or below that rank
# return a data.frame with those counts.  for example
#
# head(tf.partner.distribution(tbl.models, "TAL1", 2))
#    tf.1  tf.2 count
#  1 TAL1  TAL1    21    # 21 models in which TAL1 ranks <= 2
#  2 TAL1  IRF5     4    # there are 4 models in which TAL1 and IRF5 both rank <= 2
#  3 TAL1 THAP1     3
#  4 TAL1  JDP2     2
#  5 TAL1  ATF1     1
#  6 TAL1  FOSB     1
#
tf.partner.distribution <- function(tbl.models, tf, rank.max)
{
   tbl.tf.top <- subset(tbl.models, tf.symbol==tf & rank <= rank.max)
   tf.targets <- tbl.tf.top$target.symbol
   tbl.targets <- subset(tbl.models, target.symbol %in% tf.targets & rank <= rank.max)
   tbl.partners.dist <- as.data.frame(table(tbl.targets$tf.symbol), stringsAsFactors=FALSE)
   tbl.out <- tbl.partners.dist[order(tbl.partners.dist$Freq, decreasing=TRUE),]
   tbl.out$tf.1 <- tf
   colnames(tbl.out) <- c("tf.2", "count", "tf.1")
   rownames(tbl.out) <- NULL
   tbl.out[, c("tf.1", "tf.2", "count")]

} # tf.partner.distribution
#------------------------------------------------------------------------------------------------------------------------
test_tf.partner.distritubion <- function()
{
   printf("--- test_tf.partner.distritubtion")
   tbl.dist <- head(tf.partner.distribution(tbl.models, "TAL1", 2))
   checkEquals(tbl.dist$tf.2, c("TAL1", "IRF5", "THAP1", "JDP2", "ATF1", "FOSB"))
   checkEquals(tbl.dist$count, c(21, 4, 3, 2, 1,  1))

} # test_tf.partner.distribution
#------------------------------------------------------------------------------------------------------------------------
# get all direct interactions betweeen gene1 and gene2.  slow!

psicquic.pair <- function(gene1, gene2, quiet=TRUE)
{
   interactionType <- "direct interaction"
   tbl.i <- interactions(psicquic, id=c(gene1, gene2), species="9606", speciesExclusive=TRUE,  quiet=quiet,
                         type=interactionType)
   if(nrow(tbl.i) == 0)
      return(data.frame())

   dim(tbl.i)
   tbl.i <- addGeneInfo(psicquic.id.mapper, tbl.i)
   all.partners <- sort(unique(c(tbl.i$A.name, tbl.i$B.name)))
   subset(tbl.i, A.name %in% all.partners & B.name %in% all.partners)[,c("A.name", "B.name", "type", "detectionMethod",
                                                                         "confidenceScore", "provider", "firstAuthor")]

} # psicquic.pair
#------------------------------------------------------------------------------------------------------------------------
direct.interactors <- function(gene, quiet=TRUE)
{
   interactionType <- "direct interaction"
   tbl.i <- interactions(psicquic, id=gene, species="9606", speciesExclusive=TRUE,  quiet=quiet,
                         type=interactionType)
   if(nrow(tbl.i) == 0)
      return(c())
   tbl.i <- addGeneInfo(psicquic.id.mapper, tbl.i)
   all.partners <- sort(unique(c(tbl.i$A.name, tbl.i$B.name)))

   if(gene %in% all.partners)
      all.partners <- all.partners[-grep(gene, all.partners)]

   if("-" %in% all.partners)
      all.partners <- all.partners[-grep("-", all.partners, fixed=TRUE)]

   all.partners

} # direct.interactors
#------------------------------------------------------------------------------------------------------------------------

tfs.le.5 <- unique(subset(tbl.models, rank <= 5)$tf.symbol)  #  494
tfs.le.2 <- unique(subset(tbl.models, rank <= 2)$tf.symbol)  #  477

tfs.this.run <- tfs.le.5

tbls.all <- lapply(tfs.this.run, function(tf) tf.partner.distribution(tbl.models, tf, 5))

tbl.combined <- do.call(rbind, tbls.all)
tbl.combined$tf.1 <- as.character(tbl.combined$tf.1)
tbl.combined$tf.2 <- as.character(tbl.combined$tf.2)

#--------------------------------------------------------------------------------
# sample exploration: of the high-ranking tfs found in MEF2C-involved models
# do any have known physical interactions?  this is a weak examples
#--------------------------------------------------------------------------------
tbl.study <- subset(tbl.combined, tf.1 == "MEF2C")
possible.interactors <- tbl.study$tf.2[-1]
mef2c.known.interactors <- direct.interactors("MEF2C", quiet=TRUE)
intersect(possible.interactors, mef2c.known.interactors)   # just SMAD2 and SP1, in 7 and 19 models respectively
psicquic.pair("SMAD2", "MEF2C")  # a pull down reference, and some others

