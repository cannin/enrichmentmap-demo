#' Matrix
#' 
#' @param n entries
#' @param r number of rows 
#' @param c number of columns 
#' 
#' @concept enrichmentmap
#' @export
getMatrix <- function(n, r, c) {
  t1 <- runif(n)
  t2 <- matrix(t1, nrow=r, ncol=c)
  
  return(t2)
}

#' JSON
#' 
#' @param saveFile TBA
#' 
#' @concept enrichmentmap 
#' @export
getJson <- function(saveFile=TRUE) {
  t1 <- system.file("example.json", package="enrichmentmap")

  con <- file(t1)
  t2 <- readLines(con)
  close(con)
  
  if(saveFile) {
    fileConn <- file("example.json")
    writeLines(t2, fileConn)
    close(fileConn)    
  }

  return(t2)
}

#' JSON
#' 
#' @param p a param 
#' 
#' @concept enrichmentmap 
#' @export
getJsonStr <- function(p) {
  if(p == 0) {
    t1 <- system.file("cy.json", package="enrichmentmap")
  } else {
    t1 <- system.file("cy2.json", package="enrichmentmap")
  }
  
  con <- file(t1)
  t2 <- readLines(con)
  close(con)
  
  t3 <- paste0(t2, collapse="")
  
  return(t3)
}

#' JSON
#' 
#' @concept enrichmentmap 
#' @export
getSimple <- function() {
  return(2)
}

