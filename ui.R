library(shiny)
library(DT)
library(rcytoscapejs)

# initialGeneSet <- "BRAF\nALK\nROS1\nHRAS\nNRAS\nNF1\nFLT3\nKRAS\nEGFR\nMYC\nMET\nROS\nNTRK1\nMAP2K2\nMAP2K1\nIDH1\nBCL2\nIDH2\nBRCA1\nBRCA2\nMLL\nTP53\nRARA\nPDGFRB\nPDGFRA\nRET\nKIT"

set.seed(22)
initialGeneSet <- read.table("gene_list.txt", header=FALSE, stringsAsFactors = FALSE)
initialGeneSet <- initialGeneSet$V1
initialGeneSet <- sample(initialGeneSet, 50)
initialGeneSet <- paste(initialGeneSet, collapse="\n")

shinyUI(
  navbarPage("PC EnrichmentMap-Lite",
             theme = shinythemes::shinytheme("flatly"),
             tabPanel("Run Analysis",
                      sidebarLayout(
                        sidebarPanel(
                          width=3,
                          textAreaInput("geneList", "Gene List", initialGeneSet, 
                                        rows=16, resize='vertical'),
                          numericInput("curCluster", "Current Enrichment Map Cluster", value=1, min=1),
                          uiOutput("clusterPathways")
                          #fileInput('file1', 'Upload Custom Gene Sets',
                          #          accept=c('text/plain'))
                        ),
                        mainPanel(
                          fluidRow(
                            column(3,
                                   h4("Recurrent Terms"),
                                   tableOutput("rTable")
                            ),
                            column(3,
                                   h4(textOutput("eCount")),
                                   tableOutput("pTable")
                            ),
                            column(6,
                                   h4("Enrichment Plot"),
                                   plotOutput("ePlot", width = '400px')
                                   #img(src='enrich.png', width='400px')
                            )
                          ),
                          fluidRow(
                            column(6,
                                   h4("Enrichment Map"),
                                   rcytoscapejsOutput("eMap", height="400px")
                            ),
                            column(6,
                                   h4("Current Pathway"),
                                   rcytoscapejsOutput("path", height="400px")
                            )
                          )
                        )
                      )
             ),
             tabPanel("Run in R",
                      #includeHTML("www/files/tmp.nb.html")
                      includeMarkdown("www/files/gsea_analysis.md")
             ),
             tabPanel("About",
                      includeMarkdown("www/files/about.md")
             )
  )
)