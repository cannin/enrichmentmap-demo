library(shiny)
library(DT)
library(tibble)
library(rcytoscapejs)

library(clusterProfiler)
library(paxtoolsr)
library(tm)
library(igraph)
library(magrittr)
library(dplyr)
library(jsonlite)

# GENES 
set.seed(22)
initialGeneSet <- read.table("gene_list.txt", header=FALSE, stringsAsFactors = FALSE)
initialGeneSet <- initialGeneSet$V1
initialGeneSet <- sample(initialGeneSet, 50)
initialGeneSet <- paste(initialGeneSet, collapse="\n")

# GSA
# pc <- readGmt("PathwayCommons.8.Reactome.GSEA.hgnc.gmt", removePrefix = TRUE)
# pcGmt <- data.frame(ont=character(0), gene=character(0))
# for(i in 1:length(names(pc))) {
#   x <- names(pc)[i]
#   pcGmt <- rbind(pcGmt, data.frame(ont=x, gene=pc[[x]]))  
# }
#pcGmt <- saveRDS(pcGmt, "pcGmt.rds")
#head(pcGmt)
pcGmt <- readRDS("pcGmt.rds")

# NETWORK
source("toCytoscape.R")

# id <- c("Jerry", "Elaine", "Kramer", "George")
# name <- id
# nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
# source <- c("Jerry", "Jerry", "Jerry", "Elaine", "Elaine", "Kramer", "Kramer", "Kramer", "George")
# target <- c("Elaine", "Kramer", "George", "Jerry", "Kramer", "Jerry", "Elaine", "George", "Jerry")
# edgeData <- data.frame(source, target, stringsAsFactors=FALSE)

rSifSplit <- readRDS("reactomeSifSplit.rds")
rSifDf <- readRDS("reactomeSifDf.rds")

shinyServer(function(input, output, session) {
  getDat <- reactive({
    # DEBUG 
    #genes <- read.table("~/Dropbox/drug_target_tmp/clinical_trial_paper/gene_list.txt", stringsAsFactors = FALSE)
    #genes <- genes$V1
    
    genes <- input$geneList
    genes <- strsplit(genes, "\n")[[1]]
    cat(str(genes))
    
    # Has all need information
    egmt <- enricher(genes, TERM2GENE=pcGmt)
    
    return(egmt)
  })
  
  getEnrichMap <- reactive({
    egmt <- getDat()
    curCluster <- input$curCluster
    
    g <- enrichMap(egmt, n=10, vertex.label.font=0.05)
    ebc <- cluster_edge_betweenness(g)
    e1 <- membership(ebc)
    e2 <- names(which(e1 == curCluster))
    
    cat("Max Clusters: ", max(e1), "\n")
    
    tmp <- list(g=g, clusterCount=max(e1), clusterMembers=e2, curCluster=curCluster)
    
    return(tmp)
  })
  
  output$rTable <- renderTable({
    tmp <- getEnrichMap()
    e2 <- tmp$clusterMembers
    
    e3 <- gsub("[[:punct:]]", " ", e2)
    corpus <- Corpus(VectorSource(e3))
    corpus <- tm_map(corpus, PlainTextDocument)
    #corpus <- tm_map(corpus, removePunctuation, preserve_intra_word_dashes=TRUE)
    corpus <- tm_map(corpus, removeWords, stopwords('english'))
    #corpus <- tm_map(corpus, removeWords, c("events", "synthesis", "causes", "pathway", "superpathway", "regulation", "process"))
    dtm <- TermDocumentMatrix(corpus)
    freq <- colSums(as.matrix(dtm))
    findFreqTerms(dtm)
    
    terms <- findFreqTerms(dtm)
    termFreq <- dtm[terms,] %>%
      as.matrix() %>%
      rowSums()  %>% 
      data.frame(Term = terms, Frequency = .) %>%  
      arrange(desc(Frequency))
    dat <- termFreq[termFreq$Frequency > 1,]
    
    # dat <- tribble(
    #   ~word, ~frequency,
    #   "a", 2,
    #   "b", 1,
    #   "b", 1,
    #   "b", 1,
    #   "b", 1
    # )
    
    colnames(dat) <- c("Term", "Count")
    
    if(nrow(dat) > 10) {
      return(dat[1:10,])
    }
    
    dat
  })

  output$pTable <- renderTable({
    egmt <- getDat()
    tmp <- getEnrichMap()
    e2 <- tmp$clusterMembers
    
    idx <- which(egmt$ID %in% e2)
    z1 <- strsplit(egmt$geneID[idx], "/")
    z2 <- sort(table(unlist(z1)), decreasing = TRUE)
    dat <- data.frame(a=names(z2), b=z2, stringsAsFactors = FALSE)
    
    if(ncol(dat) == 3) {
      dat <- dat[,2:3]
    }
    
    colnames(dat) <- c("Gene", "Count")
    
    cat("pT: ", str(dat), "\n")
    
    # dat <- tribble(
    #   ~word, ~frequency,
    #   "a", 2,
    #   "b", 1
    # )

    #DT::datatable(dat, rownames=FALSE, style="bootstrap", selection="none", escape=FALSE, filter = 'top')
    
    if(nrow(dat) > 10) {
     return(dat[1:10,])
    }
    
    dat
  })
  
  output$clickedNode = renderText({
    input$clickedNode
  })
  
  output$eCount <- renderText({
    tmp <- getEnrichMap()
    a <- tmp$clusterCount
    b <- tmp$curCluster
    
    t2 <- paste0("Cluster Genes (Cluster: ", b, " Total: ", a, ")")
    t2
  })

  output$ePlot <- renderPlot({
    egmt <- getDat()
    #head(egmt)
    barplot(egmt, showCategory=10, colorBy="p.adjust", font.size=8)
  })
    
  output$eMap <- renderRcytoscapejs({
    egmt <- getDat()
    tmp <- getEnrichMap()
    g <- tmp$g
    
    id <-  V(g)$name
    name <- id
    nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
    
    edgeData <- get.edgelist(g)
    edgeData <- as.data.frame(edgeData)
    colnames(edgeData) <- c("source", "target")
    
    cyNetwork <- createCytoscapeJsNetwork(nodeData, edgeData, nodeLabelColor="#000000")
    #tmp <- fromJSON("enrich.json", simplifyVector = FALSE)
    #cyNetwork <- tmp$elements
    rcytoscapejs(nodeEntries=cyNetwork$nodes, edgeEntries=cyNetwork$edges)
  })
  
  output$path <- renderRcytoscapejs({
    #curPathway <- names(rSifSplit)[1]
    #curPathway <- "RAF activation"
    #curPathway <- input$clickedNode
    curPathway <- input$curPathway
    
    if(is.null(curPathway)) {
      curPathway <- "RAF activation"
    }
    
    #DEBUG
    #cat("CP: ", str(curPathway), "\n")
    
    idx <- rSifSplit[[curPathway]]
    curPathwayEdges <- rSifDf$edges[idx, 1:3]
    intTypes <- getSifInteractionCategories()
    curPathwayEdges <- filterSif(curPathwayEdges, interactionTypes = intTypes$BetweenProteins)
    
    cat("CP: ", str(curPathwayEdges), "\n")
    
    if(nrow(curPathwayEdges) > 300) {
      curPathwayEdges <- curPathwayEdges[1:300, ]
    }
    
    id <- unique(c(curPathwayEdges$PARTICIPANT_A, curPathwayEdges$PARTICIPANT_B))
    name <- id
    nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
    edgeData <- curPathwayEdges[, c(1,3)]
    edgeData <- unique(edgeData)
    colnames(edgeData) <- c("source", "target")
    
    cyNetwork <- createCytoscapeJsNetwork(nodeData, edgeData)
    rcytoscapejs(nodeEntries=cyNetwork$nodes, edgeEntries=cyNetwork$edges)
  })
  
  output$clusterPathways <- renderUI({
    tmp <- getEnrichMap()
    e2 <- tmp$clusterMembers
    
    tagList(
      selectInput("curPathway", "Current Cluster Pathway", e2, multiple = FALSE, selectize = TRUE)
    )
  })
})