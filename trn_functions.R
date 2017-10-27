# Examples of how to explore genome-scale TRNs

library(dplyr)
library(data.table)

# functions to be used
#---------------------------------------------------------------------------
# pull out the top TFs for a given target gene

getTarget <- function(trn, geneA)
{
	temp <- subset(trn, target.gene == geneA)
	temp[order(temp$pcaMax, decreasing=TRUE),]
}

#----------------------------------------------------------------------------
# pull out the top targets for a given TF (without rank)
getTF <- function(trn, geneA)
{
	temp <- subset(trn, gene == geneA)
	temp[order(temp$pcaMax, decreasing=TRUE),]
}

#----------------------------------------------------------------------------
# give a rank of each TF for the target gene
trnTFRank <- function(trn)
{
trn <- data.frame(trn %>% group_by(target.gene) %>% mutate(Rank = rank(-pcaMax)))
}

#----------------------------------------------------------------------------
# under the hood function
# Function for pulling out genes + rank
pullGeneAndRank <- function(my.gene,df){

    df %>% filter(target.gene == my.gene) %>%
    mutate(rank = rank(-pcaMax)) %>%
    dplyr::select(gene, rank)
}

#----------------------------------------------------------------------------
# function that takes a list of genes to see what TFs are top regulators
tfEnrich <- function(trn, gene.list)
{
df.list <- lapply(gene.list, pullGeneAndRank, trn)
all.df <- rbindlist(df.list)

# Create summary statistics
final.df <- all.df %>% group_by(gene) %>%
summarise(frequency = n(), avg.rank = mean(rank), sd.rank = sd(rank),
top.rank = min(rank), bot.rank = max(rank)) %>%
arrange(desc(frequency),avg.rank)
final.df <- as.data.frame(final.df)
arrange(head(final.df, 50), avg.rank)
}

#----------------------------------------------------------------------------
# load the data
#skin <- readRDS("primary.trn.rds")
#fibro <- readRDS("fibroblast.trn.rds")
#gene.list <- scan("ecm_list2", what="", sep="\n")

#----------------------------------------------------------------------------
# add ranks to the data
#skin <- trnTFRank(skin)
#fibro <- trnTFRank(fibro)

# pull out the top regulators of COL1A1
#getTarget(skin, "COL1A1")

# pull out the top targets of CREB3L1 (get rid of everything after the pipe (%>% if you don't want a full list)
#getTF(skin, "CREB3L1") %>% filter(Rank <= 3)

# check a list of ECM genes to see which TFs drive expression
#tfEnrich(skin, gene.list)

# get the CREB3L1 targets that are in the ecm gene list
#getTF(skin, "CREB3L1") %>% filter(target.gene %in% gene.list) %>% filter(Rank <= 3)
