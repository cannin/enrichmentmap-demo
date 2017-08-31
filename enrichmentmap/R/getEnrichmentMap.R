#' Get Enrichment Map
#' 
#' @param genes a vector of genes
#' @param curCluster integer, the current cluster 
#' @param curPathway character, the current pathway 
#' 
#' @return a list of outputs 
#' 
#' @example 
#' \dontrun{
#' tmp <- getEnrichmentMap(c("BRAF", "TP53", "MDM2")) 
#' }
#' 
#' @export 
getEnrichmentMap <- function(genes, curCluster=1, curPathway="RAF activation") {
  # DEMO ----
  #genes <- c("BRAF", "TP53", "MDM2")
  #curCluster <- 1
  #curPathway <- "RAF activation"
  
  # LOAD DATA ----
  pcGmtPath <- system.file("extdata", "pcGmt.rds", package="enrichmentmap")
  pcGmt <- readRDS(pcGmtPath)
  
  reactomeSifPath <- system.file("extdata", "reactomeSif.rds", package="enrichmentmap")
  reactomeSif <- readRDS(reactomeSifPath)
  
  reactomeSifDfPath <- system.file("extdata", "reactomeSifDf.rds", package="enrichmentmap")
  reactomeSifDf <- readRDS(reactomeSifDfPath)
  
  reactomeSifSplitPath <- system.file("extdata", "reactomeSifSplit.rds", package="enrichmentmap")
  reactomeSifSplit <- readRDS(reactomeSifSplitPath)
  
  # GET ENRICHMENT ----
  egmt <- enricher(genes, TERM2GENE=pcGmt)
  
  # GET SPECIFIC CLUSTER ----
  g <- enrichMap(egmt, n=10, vertex.label.font=0.05)
  ebc <- cluster_edge_betweenness(g)
  e1 <- membership(ebc)
  clusterMembers <- names(which(e1 == curCluster))
  
  # TMP
  cat("Max Clusters: ", max(e1), "\n")
  
  # GET TERM-COUNT ----
  # termCountDf <- tryCatch({
  #   e3 <- gsub("[[:punct:]]", " ", clusterMembers)
  #   corpus <- Corpus(VectorSource(e3))
  #   corpus <- tm_map(corpus, PlainTextDocument)
  #   #corpus <- tm_map(corpus, removePunctuation, preserve_intra_word_dashes=TRUE)
  #   corpus <- tm_map(corpus, removeWords, stopwords('english'))
  #   #corpus <- tm_map(corpus, removeWords, c("events", "synthesis", "causes", "pathway", "superpathway", "regulation", "process"))
  #   dtm <- TermDocumentMatrix(corpus)
  #   freq <- colSums(as.matrix(dtm))
  #   findFreqTerms(dtm)
  # 
  #   terms <- findFreqTerms(dtm)
  #   termFreq <- dtm[terms,] %>%
  #     as.matrix() %>%
  #     rowSums()  %>%
  #     data.frame(Term = terms, Frequency = .) %>%
  #     arrange(desc(Frequency))
  #   termCountDf <- termFreq[termFreq$Frequency > 1,]
  # 
  #   colnames(termCountDf) <- c("Term", "Count")
  # 
  #   if(nrow(termCountDf) > 10) {
  #     termCountDf <- termCountDf[1:10,]
  #   }
  # 
  #   termCountDf
  # }, error = function(e) {
  #   #cat("ERROR: ", e, "\n")
  #   message(e)
  # 
  #   termCountDf <- NULL
  # })
  
  termCountDf <- data.frame(Term="TERM", Count=0, stringsAsFactors = FALSE)
  
  # GET GENE-COUNT DF ----
  idx <- which(egmt$ID %in% clusterMembers)
  z1 <- strsplit(egmt$geneID[idx], "/")
  z2 <- sort(table(unlist(z1)), decreasing = TRUE)
  geneCountDf <- data.frame(a=names(z2), b=z2, stringsAsFactors = FALSE)
  
  if(ncol(geneCountDf) == 3) {
    geneCountDf <- geneCountDf[,2:3]
  }
  
  colnames(geneCountDf) <- c("Gene", "Count")
  
  cat("pT: ", str(geneCountDf), "\n")

  if(nrow(geneCountDf) > 10) {
    geneCountDf <- geneCountDf[1:10,]
  }
  
  # GET ENRICHMENTMAP ----
  id <-  V(g)$name
  name <- id
  nodeData <- data.frame(id, name, stringsAsFactors=FALSE)
  
  edgeData <- get.edgelist(g)
  edgeData <- as.data.frame(edgeData)
  colnames(edgeData) <- c("source", "target")
  
  ## Fix this??? 
  eMapCyNetwork <- createCytoscapeJsNetwork(nodeData, edgeData, nodeLabelColor="#000000")
  #tmp <- fromJSON("enrich.json", simplifyVector = FALSE)
  #cyNetwork <- tmp$elements
  #rcytoscapejs(nodeEntries=cyNetwork$nodes, edgeEntries=cyNetwork$edges)
  
  # GET CURRENT PATHWAY ----
  idx <- reactomeSifSplit[[curPathway]]
  tmpCurPathwayEdges <- reactomeSifDf$edges[idx, 1:3]
  #intTypes <- getSifInteractionCategories()
  #curPathwayEdges <- filterSif(curPathwayEdges, interactionTypes = intTypes$BetweenProteins)
  
  intTypes <- c("controls-state-change-of", "controls-expression-of", "controls-degradation-of", 
    "controls-transport-of", "catalysis-precedes", "in-complex-with")
  idx <- which(tmpCurPathwayEdges$INTERACTION_TYPE %in% intTypes)
  curPathwayEdges <- tmpCurPathwayEdges[idx,]
  
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
  
  curPathwayCyNetwork <- createCytoscapeJsNetwork(nodeData, edgeData)
  
  # CREATE RESULTS ----
  tmpEgmtDf <- as.data.frame(egmt)
  results <- list(g=g, egmt=tmpEgmtDf, clusterCount=max(e1), clusterMembers=clusterMembers, curCluster=curCluster, 
                  termCountDf=termCountDf, geneCountDf=geneCountDf, eMapCyNetwork=eMapCyNetwork, 
                  curPathwayCyNetwork=curPathwayCyNetwork)  

  # SAVE TO JSON ----
  ## No asJSON for igraph
  results$g <- NULL
  
  tmp <- toJSON(results, pretty = TRUE)
  fileConn <- file("tmp.json")
  writeLines(tmp, fileConn)
  close(fileConn)
  
  return(results)
}
