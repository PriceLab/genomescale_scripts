# identify the top targets for a given TF, and give the TF rank relative to other TFs for the target gene.

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
	subset(trn, target.gene == geneA)
}

#----------------------------------------------------------------------------
# pull out the top targets for a given TF (without rank)
getTF <- function(trn, geneA)
{
	subset(trn, gene == geneA)
}
