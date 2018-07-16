# symToGeneID.R:  use org.Hs.egSYMBOL2EG and org.Hs.egALIAS2EG
# to map gene symbols to geneIDs
#------------------------------------------------------------------------------------------------------------------------
library(org.Hs.eg.db)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
assignGeneIDs <- function(symbols)
{
    geneIDs <- mget(symbols, org.Hs.egSYMBOL2EG, ifnotfound=NA)

    if(all(!is.na(geneIDs)) & all(sapply(geneIDs, length) == 1))
        return(list(mapped=geneIDs, failures=NULL, multiples=NULL))

    multiples <- names(which(sapply(geneIDs, length) > 1))

    if(length(multiples) == 0)
        multiples <- NULL

    unmapped <- names(which(is.na(geneIDs)))


    if(length(unmapped) == 0){
        unmapped <- NULL
        }
    else{
        aliasIDs <- mget(unmapped, org.Hs.egALIAS2EG, ifnotfound=NA)
        unmapped <- names(which(is.na(aliasIDs)))
        #browser("which")
        new.multiples <- names(which(sapply(aliasIDs, length) > 1))
        if(length(new.multiples) > 0)
            multiples <- c(multiples, new.multiples)
        #browser("assign")
        geneIDs[names(aliasIDs)] <- aliasIDs
        }

    if(length(unmapped) == 0)
        unmapped <- NULL
    if(all(is.na(geneIDs)))
        geneIDs <- NULL
    if(length(multiples) == 0)
        multiples <- NULL

    return(list(mapped=geneIDs, failures=unmapped, multiples=multiples))

} # assignGeneIDs
#------------------------------------------------------------------------------------------------------------------------
test_assignGeneIDs <- function()
{
    print("--- test_assignGeneIDs")

    syms.mappedInSYMBOL2EG <- c("AAK1", "AATK", "ABCA3")
    syms.mappedOnlyInALIAS2EG <- c("TBCKL", "PCTK3")
    syms.multipleInSYMBOL2EG <- c("TEC", "HBD")
    syms.multipleInALIAS2EG <- c("RAGE", "IL8RA")
    syms.unknown <- c("SgK269", "bogus-gene-name", "Sep-7")


    x <- assignGeneIDs(syms.mappedInSYMBOL2EG)
    checkEquals(names(x$mapped), syms.mappedInSYMBOL2EG)
    checkEquals(as.character(x$mapped), c("22848", "9625", "21"))
    checkTrue(is.null(x$failures))
    checkTrue(is.null(x$multiples))

    x <- assignGeneIDs(syms.unknown)
    checkTrue(is.null(x$mapped))
    checkEquals(x$failures, syms.unknown)
    checkTrue(is.null(x$multiples))

    x <- assignGeneIDs(syms.multipleInSYMBOL2EG)
    checkEquals(names(x$mapped), syms.multipleInSYMBOL2EG)
    checkTrue(all(sapply(x$mapped, length) > 1))
    checkTrue(is.null(x$failures))
    checkEquals(x$multiples, syms.multipleInSYMBOL2EG)

    #x <- assignGeneIDs(syms)
    #checkEquals(x$failures, "SgK269")
    #checkEquals(x$multiples, "TEC")

    x <- assignGeneIDs(syms.mappedOnlyInALIAS2EG)
    checkEquals(names(x$mapped), syms.mappedOnlyInALIAS2EG)
    checkEquals(unlist(x$mapped, use.names=FALSE), c("93627", "5129"))
    checkTrue(is.null(x$failures))
    checkTrue(is.null(x$multiples))

       # now try some combinations.  first, 2 mapped in alias, 3 completely unknown, 2 multiple
    syms <- c(syms.mappedOnlyInALIAS2EG, syms.unknown, syms.multipleInSYMBOL2EG)
    x <- assignGeneIDs(syms)
    checkEquals(names(x$mapped), syms)
    checkEquals(x$failures, syms.unknown)
    checkEquals(x$multiples, syms.multipleInSYMBOL2EG)

    syms <- c(syms.mappedInSYMBOL2EG,
              syms.mappedOnlyInALIAS2EG,
              syms.multipleInSYMBOL2EG,
              syms.multipleInALIAS2EG,
              syms.unknown)

    x <- assignGeneIDs(syms)
    checkEquals(names(x$mapped), syms)
    checkEquals(x$failures, syms.unknown)

    checkEquals(sort(x$multiples), sort(c(syms.multipleInSYMBOL2EG,
                                   syms.multipleInALIAS2EG)))
       # all symbols reported in multiples should indeed have
       # multiple entries in mapped

    checkTrue(all(lapply(x$mapped[x$multiples], length) > 1))


} # test_assignGeneIDs
#------------------------------------------------------------------------------------------------------------------------
