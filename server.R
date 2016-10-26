
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyBS)
library(DT)
library(fgsea)
library(ggplot2)

pathways <- NULL
ranks <- NULL
pathway <- NULL
imagename <- NULL

update_plot <- function() {
    plot <- plotEnrichment(pathways[[pathway]], ranks) + labs(title=pathway)
    imagename <<- paste0('www/', pathway, '_plot.png')
    ggsave(filename = imagename, plot = plot)
}

get_pathways <- function(gmt.file) {
    pathwayLines <- strsplit(readLines(gmt.file), "\t")
    aux <- lapply(pathwayLines, tail, -2)
    names(aux) <- sapply(pathwayLines, head, 1)
    return(aux)
}

get_ranks <- function(rnk.file) {
    aux <- read.table(rnk.file, header=T, colClasses = c("character", "numeric"))
    aux <- setNames(aux$t, aux$ID)
    return(aux)
}

createLink <- function(val) {
    sprintf('<button id="%s" onclick="btnclick(this.id)">%s</button>', val, 'Show plot')
}

shinyServer(function(input, output) {

    output$downloadExampleRNK <- downloadHandler(
        filename = 'sample.rnk',
        content = function(con) {
            file.copy('./data/sample.rnk', con)
        }
    )

    output$downloadExampleGMT <- downloadHandler(
        filename = 'sample.gmt',
        content = function(con) {
            file.copy('./data/sample.gmt', con)
        }
    )

    output$contents <- renderDataTable(
        {
            rnk.file <- input$rnkfile$datapath
            gmt.file <- input$gmtfile$datapath

            if (is.null(rnk.file) | is.null(gmt.file))
                return(NULL)

            # hide sidebar and show loader animation
            js$hideSidebar()
            shinyjs::show("container")

            # load and process data
            ranks <<- get_ranks(rnk.file)
            pathways <<- get_pathways(gmt.file)
            res <- NULL
            res <- fgsea(pathways, ranks, nperm=10000, maxSize=500)
            res$plot <- createLink(res$pathway)
            res$pval <- format(round(res$pval, 6), nsmall = 6)
            res$padj <- format(round(res$padj, 6), nsmall = 6)
            res$ES <- format(round(res$ES, 6), nsmall = 6)
            res$NES <- format(round(res$NES, 6), nsmall = 6)
            res[,leadingEdge:=NULL]
            sapply(res$pathway, function(x) {onclick(x, update_plot(x))})

            # hide loader and show table
            shinyjs::hide("container")
            res
        },
        rownames = FALSE,
        escape = FALSE,
        selection = 'none',
        options = list(
            lengthChange = TRUE,
            pageLength = 20,
            lengthMenu = c(10, 20, 50, 100, -1),
            order = list(list(5, 'desc'), list(4, 'desc'))
        )
    )

    observe({
        new_pathway <- input$JsMessage
        if (!is.null(new_pathway)) {
            if (is.null(pathway)) {
                pathway <<- new_pathway
                update_plot()
            }
            if (new_pathway != pathway) {
                pathway <<- new_pathway
                update_plot()
            }
        }
    })
})
