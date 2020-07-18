library(shiny)
library(shinyModules)
library(EnsDb.Hsapiens.v79)
library(TrenaProjectHG38.generic)
library(GenomicRanges)
#----------------------------------------------------------------------------------------------------
load("tbl.twE12.deg.RData")
ensg.deg.pos <- subset(tbl.deg, log2FoldChange >  8 & -log10(pvalue) > 8)$Gene_ID
sym.deg.pos <- select(EnsDb.Hsapiens.v79, key=ensg.deg.pos,  columns=c("SYMBOL", "GENEID"),
                      keytype="GENEID")$SYMBOL
printf <- function(...) print(noquote(sprintf(...)))
printf("somewhat arbitrary selection of %d deg + genes", length(sym.deg.pos))

#----------------------------------------------------------------------------------------------------
ui <- fluidPage(

  titlePanel("Explore ChIP-seq Binding of Some DEG+ Genes"),
  sidebarLayout(
    sidebarPanel(
      selectInput("degPositiveGenesSelector",
                  "Choose Gene of Interest",
                  c(" - ", head(sym.deg.pos))),
      width=2
    ),

    mainPanel(
       div(igvUI("igv"),
           style="height: 2000px; margin: 10px; margin-bottom: 5px; padding: 10px;"),
       style="margin-top:5px;",
       width=10
       )
  )

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
      if(selectedGene != " - ")
         showGenomicRegion(session, selectedGene)
      })


} # server
#----------------------------------------------------------------------------------------------------
runApp(shinyApp(ui, server), host="0.0.0.0", port=7890)

