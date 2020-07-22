library(shiny)
library(shinyModules)
library(EnsDb.Hsapiens.v79)
library(ghdb)

#tp <- TrenaProjectHG38.generic()
#library(TrenaProjectHG38.generic)

library(GenomicRanges)
ghdb <- GeneHancerDB()

#----------------------------------------------------------------------------------------------------
load("tbl.twE12.deg.RData")
ensg.deg.pos <- subset(tbl.deg, log2FoldChange >  8 & -log10(pvalue) > 8)$Gene_ID
sym.deg.pos <- select(EnsDb.Hsapiens.v79, key=ensg.deg.pos,  columns=c("SYMBOL", "GENEID"),
                      keytype="GENEID")$SYMBOL
printf <- function(...) print(noquote(sprintf(...)))
printf("somewhat arbitrary selection of %d deg + genes", length(sym.deg.pos))

tbl.peaks.twe12 <- read.table("twe12.bed", sep="\t", skip=1, header=FALSE, as.is=TRUE)[, c(1,2,3,5)]
tbl.peaks.twtw <- read.table("twtw.bed", sep="\t", skip=1, header=FALSE, as.is=TRUE)[, c(1,2,3,5)]
colnames(tbl.peaks.twe12) <- c("chrom", "start", "end", "score")
colnames(tbl.peaks.twtw) <- c("chrom", "start", "end", "score")
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
                               width=2
                           ),
                           mainPanel(
                               div(igvUI("igv"),
                                   style="height: 2000px; margin: 10px; margin-bottom: 5px; padding: 10px;"),
                               style="margin-top:5px;",
                               width=10
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
             tbl.sub.twe12 <- subset(tbl.peaks.twe12, chrom==loc.chrom & start > loc.start & end < loc.end)
             tbl.sub.twtw <- subset(tbl.peaks.twtw, chrom==loc.chrom & start > loc.start & end < loc.end)
             loadBedGraphTrack(session, "TWE12", tbl.sub.twe12, autoscale=TRUE, color="random")
             loadBedGraphTrack(session, "TWTW", tbl.sub.twtw, autoscale=TRUE, color="random")

             gr.gh <- GRanges(tbl.gh)
             gr.twe12 <- GRanges(tbl.sub.twe12)
             gr.twtw <- GRanges(tbl.sub.twtw)
             tbl.ov.twe12 <- as.data.frame(findOverlaps(gr.twe12, gr.gh))
             tbl.ov.twtw <- as.data.frame(findOverlaps(gr.twtw, gr.gh))

             tbl.twe12.gh <- tbl.sub.twe12[unique(tbl.ov.twe12[, 1]),]
             loadBedGraphTrack(session, "TWE12.gh", tbl.twe12.gh, autoscale=TRUE, color="random")

             tbl.twtw.gh <- tbl.sub.twtw[unique(tbl.ov.twtw[, 1]),]
             loadBedGraphTrack(session, "TWTW.gh", tbl.twtw.gh, autoscale=TRUE, color="random")
             }
         } # if legit selectedGene
      })


} # server
#----------------------------------------------------------------------------------------------------
runApp(shinyApp(ui, server), host="0.0.0.0", port=3838)

