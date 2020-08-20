library(shiny)
library(shinyModules)
library(EnsDb.Hsapiens.v79)
library(ghdb)
library(rtracklayer)

#tp <- TrenaProjectHG38.generic()
#library(TrenaProjectHG38.generic)

library(GenomicRanges)
ghdb <- GeneHancerDB()

#----------------------------------------------------------------------------------------------------
load("tbl.twE12.deg.RData")

chip.dir <- "~/github/genomescale_scripts/bobRostomily/AM.HMRI.TWIST1.CSeq.37738.2"

ensg.deg.pos <- subset(tbl.deg, log2FoldChange >  8 & -log10(pvalue) > 8)$Gene_ID
sym.deg.pos <- select(EnsDb.Hsapiens.v79, key=ensg.deg.pos,  columns=c("SYMBOL", "GENEID"),
                      keytype="GENEID")$SYMBOL
printf <- function(...) print(noquote(sprintf(...)))
printf("somewhat arbitrary selection of %d deg + genes", length(sym.deg.pos))

tbl.peaks.twe12 <- read.table("twe12.bed", sep="\t", skip=1, header=FALSE, as.is=TRUE)[, c(1,2,3,5)]
tbl.peaks.twtw <- read.table("twtw.bed", sep="\t", skip=1, header=FALSE, as.is=TRUE)[, c(1,2,3,5)]
colnames(tbl.peaks.twe12) <- c("chrom", "start", "end", "score")
colnames(tbl.peaks.twtw) <- c("chrom", "start", "end", "score")

tbls.deg <- get(load("degGenes.RData"))

experiments <- list("dTWB98vsCNTRL98_DEG_all.tsv"="TW.del",
                    "TEvsCNTRL.T98G_DEG_all.tsv"="TW.E12",
                    "TS68AEvsCNTRL.T98G_DEG_all.tsv"="TW.mutant.vs.WT",
                    "TS68AEvsTE.T98G.DEG_all.tsv"="TW.mutant.vs.TWE12",
                    "TWmonomer.vsCNTRLT98G_DEG_all.tsv"="TW.monomer",
                    "TwTwvsCNTRL.T98G_DEG_all.tsv"="TW.dimer")
state <- new.env(parent=emptyenv())
#----------------------------------------------------------------------------------------------------
ui <- fluidPage(

  titlePanel("Explore ChIP-seq Binding of Some DEG+ Genes"),


  tabsetPanel(type = "tabs",
              tabPanel("Introduction", includeHTML("intro.html")),
              tabPanel("ChiP-SEQ in Genome Context",
                       sidebarLayout(
                           sidebarPanel(
                               selectInput("degPositiveGenesSelector",
                                           "Choose Gene of Interest",
                                           c(" - ", head(sym.deg.pos, n=20))),
                               selectInput("experimentSelector",
                                           "Experiment:",
                                           c(" - ", as.character(experiments))),
                               width=3
                           ),
                           mainPanel(
                               div(igvUI("igv"),
                                   style="height: 2000px; margin: 10px; margin-bottom: 5px; padding: 10px;"),
                               style="margin-top:5px;",
                               width=9
                           ) # mainPanel
                       ) # sidebarLayout
                       ) # chip-seq tab
              ) # tabsetPanel
) # ui
#----------------------------------------------------------------------------------------------------
server <- function(session, input, output) {

  callModule(igvServer, "igv",
             genome="hg38",
             geneModelDisplayMode="COLLAPSED",
             locus="all")
             #locus="chr19:44,854,808-44,940,011") #"APOE")

   observeEvent(input$degPositiveGenesSelector, {
      selectedGene <- input$degPositiveGenesSelector;
      if(selectedGene != " - "){
         tbl.gh <- getEnhancers(ghdb, selectedGene, maxSize=40000)
         state$gh <- tbl.gh
         if(nrow(tbl.gh) == 0){
           showNotification("No genehancer information for this gene.", type="error")
         } else {
             loc.chrom <- tbl.gh$chrom[1]
             shoulder <- 2000
             loc.start <- min(tbl.gh$start) - shoulder
             loc.end   <- max(tbl.gh$end) + shoulder
             chromLoc <- sprintf("%s:%d-%d", loc.chrom, loc.start, loc.end)
             showGenomicRegion(session, chromLoc)
             tbl.gh <- subset(tbl.gh, elite)
             loadBedGraphTrack(session, "GH", tbl.gh[, c("chrom", "start", "end", "combinedscore")],
                               autoscale=FALSE, min=0, max=50, color="brown")
         } # else
         } # if legit selectedGene
      })


    observeEvent(input$experimentSelector, {
        experiment.name <- input$experimentSelector
        if(experiment.name != " - "){
           printf("experiment: %s", experiment.name)
           chrom.loc <- state$gh$chrom[1]
           start.loc <- min(state$gh$start) - 1000
           end.loc   <- max(state$gh$end) + 1000
           gr.region <- GRanges(seqnames=chrom.loc, IRanges(start=start.loc, end=end.loc))
           success <- FALSE
           switch(experiment.name,
                "TW.dimer" = {
                    tbl.peaks <- subset(tbl.peaks.twtw, chrom==chrom.loc & start > start.loc & end < end.loc)
                    file <- file.path(chip.dir, "3_079N_00YMHMRI_T98G-TW-TW_TWIST1_hg38_i90_dmnorm_signal.bw")
                    gr.chip <- import(con=file, which=gr.region)
                    tbl.deg <- tbls.deg[["TWTW"]]
                    success <- TRUE
                    }
                 ) # switch on experiment.name
           if(success){
              tbl.chip <- as.data.frame(gr.chip)[, c("seqnames", "start", "end", "score")]
              colnames(tbl.chip)[1] <-"chrom"
              tbl.chip$chrom <- as.character(tbl.chip$chrom)
              loadBedGraphTrack(session, sprintf("%s-ChIP", experiment.name), tbl.chip,
                                autoscale=TRUE, color="red")
              loadBedGraphTrack(session, sprintf("%s-peaks", experiment.name), tbl.peaks,
                                autoscale=TRUE, color="gray")
              selectedGene <- isolate(input$degPositiveGenesSelector)
              if(selectedGene %in% tbl.deg$GeneName){
                 tbl.expression <- subset(tbl.deg, GeneName==selectedGene)
                 tbl.track <- tbl.expression[, c("Chrom", "Start", "End", "log2FoldChange")]
                 colnames(tbl.track) <- c("chrom", "start", "end", "score")
                 loadBedGraphTrack(session, sprintf("%s-RNA", experiment.name), tbl.track,
                                   autoscale=FALSE, min=-5, max=5, color="brown")
                 } # if selectedGene
              } # if success

        # tbl.sub.twe12 <- subset(tbl.peaks.twe12, chrom==loc.chrom & start > loc.start & end < loc.end)
        # tbl.sub.twtw <- subset(tbl.peaks.twtw, chrom==loc.chrom & start > loc.start & end < loc.end)
        # loadBedGraphTrack(session, "TWE12", tbl.sub.twe12, autoscale=TRUE, color="random")
        # loadBedGraphTrack(session, "TWTW", tbl.sub.twtw, autoscale=TRUE, color="random")

        # tbls.deg <- lapply(tbls, function(tbl) subset(tbl, GeneName==selectedGene))
        # track.names <- names(tbls)[unlist(lapply(tbls.deg, function(tbl.deg) nrow(tbl.deg) > 0), use.names=FALSE)]
        # tbl.deg <- do.call(rbind, tbls.deg)
        # if(nrow(tbl.deg) > 0){
        #   chrom.gene <- tbl.deg$Chrom[1]
        #   chrom.start <- tbl.deg$Start[1]
        #   chrom.end <- tbl.deg$End[1]
        #   for(r in seq_len(nrow(tbl.deg))){
        #     track.name <- sprintf("%s-RNA", track.names[r])
        #     track.score <- tbl.deg$log2FoldChange[r]
        #     tbl.graph <- data.frame(chrom=chrom.gene, start=chrom.start, end=chrom.end, score=track.score, stringsAsFactors=FALSE)
        #     track.color <- "red"
        #     if(track.score < 0) track.color <- "green"
        #     loadBedGraphTrack(session, track.name, tbl.graph, autoscale=FALSE, min=-5, max=5, color=track.color)
        #     } # for r
        #   } # if nrow(tbl.deg)
          }
        })

} # server
#----------------------------------------------------------------------------------------------------
runApp(shinyApp(ui, server), host="0.0.0.0", port=3838)
#shinyApp(ui, server)

