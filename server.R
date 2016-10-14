
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(fgsea)

get_pathways <- function(gmt.file) {
  pathwayLines <- strsplit(readLines(gmt.file), "\t")
  pathways <- lapply(pathwayLines, tail, -2)
  names(pathways) <- sapply(pathwayLines, head, 1)
  return(pathways)
}

get_ranks <- function(rnk.file) {
  ranks <- read.table(rnk.file, header=T, colClasses = c("character", "numeric"))
  ranks <- setNames(ranks$t, ranks$ID)
  return(ranks)
}

shinyServer(function(input, output) {

  output$contents <- renderTable({
    rnk.file <- input$rnkfile$datapath
    gmt.file <- input$gmtfile$datapath

    if (is.null(rnk.file) | is.null(gmt.file))
      return(NULL)

    ranks <- get_ranks(rnk.file)
    pathways <- get_pathways(gmt.file)
    res <- fgsea(pathways, ranks, nperm=10000, maxSize=500)
    res[,leadingEdge:=NULL]
    res
  })

})
