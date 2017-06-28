# identify the top targets for a given TF, and give the TF rank relative to other TFs for the target gene.
library(dplyr)
library(data.table)


addTFRank <- function(df, tf){

extractRank <- function(df,my.target,tf){

    blah <- df %>% filter(target.gene==my.target) %>% mutate(rank = rank(-pcaMax)) %>% filter(gene==tf) %>% select(rank)
    blah$rank[[1]]
}

temp1 <- filter(df, gene==tf)
temp1$rank <- NA
for(i in 1:length(temp1$gene)){

    temp1$rank[[i]] <- extractRank(df,
                                      temp1$target.gene[[i]],
                                      temp1$gene[[i]])
}

return(arrange(temp1, rank, desc(pcaMax)))

}

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

# Function for pulling out genes + rank
pullGeneAndRank <- function(my.gene,df){

    df %>% filter(target.gene == my.gene) %>%
    mutate(rank = rank(-pcaMax)) %>%
    select(gene, rank)
}

# example of applying the function
#gene.list <- scan("dev.genes", what="", sep="\n")

#df.list <- lapply(gene.list, pullGeneAndRank, mtx.all)

#all.df <- rbindlist(df.list)

# Create summary statistics
#final.df <- all.df %>% group_by(gene) %>%
#summarise(frequency = n(), avg.rank = mean(rank), sd.rank = sd(rank),
#top.rank = min(rank), bot.rank = max(rank)) %>%
#arrange(desc(frequency),avg.rank)
#final.df <- as.data.frame(final.df)
#head(final.df, 30)
