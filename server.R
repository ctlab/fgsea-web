library(data.table)
library(shiny)
library(shinydashboard)
library(shinyjs)
library(DT)
library(fgsea)
library(ggplot2)
library(svglite)
library(hash)
library(AnnotationDbi)

source("./functions.R")

if (!file.exists("data/org.annos.rda")) {
    source("prepare.R")
}
load("data/org.annos.rda")
load("data/pathways.msigdb.rda")

stripped_to_full <- NULL
full_to_stripped <- NULL

createLink <- function(val) {
    sprintf('<button id="%s" onclick="btnclick(this.id)">%s</button>', val, 'Show plot')
}

stripSpecialChars <- function(string) {
    gsub("[\\'\\;\\&\\:\\,\\.]", "", string)
}


shinyServer(function(input, output, session) {

    output$downloadExampleTSV <- downloadHandler(
        filename = 'Ctrl.vs.MandLPSandIFNg.gene.de.tsv',
        content = function(con) {
            file.copy('./data/Ctrl.vs.MandLPSandIFNg.gene.de.tsv', con)
        }
    )

    geneDEInput <- reactive({
        tsv.file <- input$tsvfile$datapath
        if (is.null(tsv.file))
            return(NULL)
        geneDe <- fread(tsv.file)
    })
    
    geneDEMetaInput <- reactive({
        geneDE <- geneDEInput()
        if (is.null(geneDE)) {
            return(NULL)
        }
        getGeneDEMeta(geneDE, org.annos)
    })

    geneRanksInput <- reactive({
        geneDEMeta <- geneDEMetaInput()
        if (is.null(geneDEMeta)) {
            return(NULL)
        }
        geneDE <- geneDEInput()
        geneRanks <- getRanks(geneDE, geneDEMeta, org.annos)
    })
    
    output$uiDEInfo <- renderUI({
        geneDEMeta <- geneDEMetaInput()
        if (is.null(geneDEMeta)) {
            return(NULL)
        }
        geneRanks <- geneRanksInput()
        
        geneRanksSym <- paste0(org.annos[[geneDEMeta$organism]]$genes[names(geneRanks), symbol], 
                               ": ", geneRanks)
        
        div(p("Detected info:"),
            p(sprintf("species: %s", geneDEMeta$organism)),
            p(sprintf("ID column: %s", geneDEMeta$columns$ID)),
            p(sprintf("ID type: %s", geneDEMeta$idType)),
            p(sprintf("Mean expression column: %s", geneDEMeta$columns$baseMean)),
            p(sprintf("Ranking column: %s", geneDEMeta$columns$stat)),
            p(sprintf("Top up genes: %s", paste0(head(geneRanksSym), collapse=", "))),
            p(sprintf("Top down genes: %s", paste0(tail(geneRanksSym), collapse=", ")))
            )
    })

    output$selectPathways <- renderUI({
        geneDEMeta <- geneDEMetaInput()
        if (is.null(geneDEMeta)) {
            return(NULL)
        }
        
        if (geneDEMeta$organism == 'Mus musculus') {
            includeHTML('static/mm-checkboxes.html')
        } else if (geneDEMeta$organism == 'Homo sapiens') {
            includeHTML('static/hs-checkboxes.html')
        } else {
            includeHTML("Unsupported organism!")
        }
    })
    
    pathwaysInput <- reactive({
        geneDEMeta <- geneDEMetaInput()
        if (is.null(geneDEMeta)) {
            return(NULL)
        }
        selectedGenesets <-input$selectedGenesets 
        if (length(selectedGenesets) == 0)
            return(NULL)
        
        pathways <- do.call(c, pathways.msigdb[[geneDEMeta$organism]][selectedGenesets])
        pathways
    })

    processButtonClick <- eventReactive(input$submitButton, {
        .messagef("Submit button clicked")
        geneDEMeta <- geneDEMetaInput()
        if (is.null(geneDEMeta)) {
            return(NULL)
        }
        
        selectedGenesets <-input$selectedGenesets 
        if (length(selectedGenesets) == 0) {
            return(NULL)
        }
        
        print(paste(selectedGenesets))
        
        geneDE <- geneDEInput()
        ranks <- geneRanksInput()
        pathways <- pathwaysInput()
        
        
        
        addClass(selector = "body", class = "sidebar-collapse")
        shinyjs::show("container")
        
        pathways <- do.call(c, pathways.msigdb[[geneDEMeta$organism]][selectedGenesets])
        
        res <- NULL
        set.seed(42)
        res <- fgsea(pathways, ranks, nperm=25000, minSize=15, maxSize=500)
        strippedPathways <- stripSpecialChars(res$pathway)
        full_to_stripped <<- hash(keys=res$pathway, values=strippedPathways)
        stripped_to_full <<- hash(keys=strippedPathways, values=res$pathway)
        res$pathway <- strippedPathways
        res$plot <- createLink(res$pathway)
        res$pval <- round(res$pval, 6)
        res$padj <- round(res$padj, 6)
        res$ES <- round(res$ES, 6)
        res$NES <- round(res$NES, 6)
        res[,leadingEdge:=NULL]
        shinyjs::hide("container")
        return(res)
    })

    output$contents <- renderDataTable(
        {
            processButtonClick()
        },
        rownames = FALSE,
        escape = FALSE,
        selection = 'none',
        options = list(
            lengthChange = TRUE,
            pageLength = 20,
            lengthMenu = c(10, 20, 50, 100, -1),
            order = list(list(2, 'asc'))
        )
    )

    observe({
        if (!is.null(input$pathway)) {
            pathways <- isolate(pathwaysInput())
            ranks <- isolate(geneRanksInput())
            pathway <- stripped_to_full[[input$pathway]]
            dir.create("www/img", showWarnings=FALSE)
            filePath <- tempfile(tmpdir = "www/img", fileext = '.svg')
            if (!file.exists(filePath)) {
                plot <- plotEnrichment(pathways[[pathway]], ranks)
                ggsave(filename = filePath, plot = plot, width=6, height=4)
            }
            # TODO: fix this and save images in tmp folder
            message = list(link = substr(filePath, 5, nchar(filePath)), pathway = full_to_stripped[[pathway]])
            print(message)
            session$sendCustomMessage(type = "imageReady", message)
        }
    })
})

