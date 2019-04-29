library(RUnit)
library(gplots)
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
   tbl.tf.top <- subset(tbl.models, tf.hgnc==tf & rank <= rank.max)
   tf.targets <- tbl.tf.top$targetGene.hgnc
   tbl.targets <- subset(tbl.models, targetGene.hgnc %in% tf.targets & rank <= rank.max)
   tbl.partners.dist <- as.data.frame(table(tbl.targets$tf.hgnc), stringsAsFactors=FALSE)
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
   checkTrue(all(c("TAL1", "IRF5",  "LYL1") %in% tbl.dist$tf.2))
   checkEquals(tbl.dist$count[1:3], c(22, 5, 2))

} # test_tf.partner.distribution
#------------------------------------------------------------------------------------------------------------------------
# collect about ten tfs, then

seed.tf <- "TAL1"
min.rank <- 5
partner.tf.must.have.count.of.at.least <- 5
goi <- subset(tf.partner.distribution(tbl.models, seed.tf, min.rank), count >= partner.tf.must.have.count.of.at.least)$tf.2
length(goi)

tbls <- lapply(goi, function(tf) tf.partner.distribution(tbl.models, tf, min.rank))
length(tbls)
tbl <- do.call(rbind, tbls)
tbl <- subset(tbl, count > 12)
dim(tbl)

tfs.all <- sort(unique(c(tbl$tf.1, tbl$tf.2)))
tf.count <- length(tfs.all)  # 192
tf.count

mtx <- matrix(0, nrow=tf.count, ncol=tf.count, dimnames=list(tfs.all, tfs.all))

for(i in seq_len(tf.count)){
   tf.1 <- tbl[i, "tf.1"]
   tf.2 <- tbl[i, "tf.2"]
   count <- tbl[i, "count"]
   mtx[tf.1, tf.2] <- count
   mtx[tf.1, tf.1] <- count
   }

diag(mtx) <- 0
heatmap.2(mtx, trace='none', col=c("#000000", rev(heat.colors(20))))
