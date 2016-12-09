
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinydashboard)
library(shinyjs)
library(DT)
library(fgsea)
library(ggplot2)
library(svglite)
library(hash)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationDbi)

pathways <- NULL
ranks <- NULL
stripped_to_full <- NULL
full_to_stripped <- NULL
detectedSpecies <- NULL
detectedFormat <- NULL
alreadyConverted <- FALSE

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

stripSpecialChars <- function(string) {
    gsub("[\\'\\;\\&\\:\\,\\.]", "", string)
}

loadPathwayFile <- function(path) {
    bonus_pathways <- get_pathways(path)
    pathways <<- c(pathways, bonus_pathways)
}

convertToEntrez <- function() {
    if (alreadyConverted | detectedFormat == 'ENTREZID') {
        return()
    }

    if (detectedSpecies == 'mm') {
        converted <- AnnotationDbi::mapIds(org.Mm.eg.db, keys=names(ranks), column="ENTREZID", keytype=detectedFormat, multiVals="first")
    } else {
        converted <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=names(ranks), column="ENTREZID", keytype=detectedFormat, multiVals="first")
    }

    print(paste('Failed to convert', sum(is.na(converted)), '/', length(converted), 'genes'))
    inds <- which(!is.na(converted))
    ranks <<- ranks[inds]
    names(ranks) <<- converted[inds]
    alreadyConverted <<- TRUE
}

detectSpecies <- function() {
    idExample <- names(ranks)[1]
    detectedFormat <<- 'ENTREZID'
    if (is.na(as.numeric(idExample))) {
        if (startsWith(idExample, 'ENS')) {
            detectedFormat <<- 'ENSEMBL'
        } else {
            detectedFormat <<- 'SYMBOL'
        }
    }
    detectedSpecies <<- tryCatch({
        converted <- AnnotationDbi::mapIds(org.Mm.eg.db, keys=idExample, column="ENTREZID", keytype=detectedFormat, multiVals="first")
        'mm'
    }, error = function(e) {
        tryCatch({
            converted <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=idExample, column="ENTREZID", keytype=detectedFormat, multiVals="first")
            'hs'
        }, error = function(e) {
            'unknown'
        })
    })
    print(paste(detectedFormat, detectedSpecies))
    return(detectedSpecies)
}

shinyServer(function(input, output, session) {

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

    output$useOwnPathwaysRadio <- renderUI({
        rnk.file <- input$rnkfile$datapath
        if (is.null(rnk.file))
            return(NULL)

        ranks <<- get_ranks(rnk.file)
        detectSpecies()
        if (detectedSpecies == 'unknown') {
            return(NULL)
        } else {
            title <- ifelse(detectedSpecies == 'hs', 'Detected species: Homo sapiens', 'Detected species: Mus musculus')
            radioButtons("useOwnPathways", title,
                         c("Use my own pathways" = TRUE,
                           "Use standard pathways" = FALSE),
                         selected = FALSE)
        }
    })

    output$selectPathways <- renderUI({

        if (is.null(input$rnkfile$datapath))
            return(NULL)

        if (is.null(detectedSpecies))
            return(NULL)

        if (detectedSpecies == 'unknown') {
            fileInput('gmtfile', 'Unknown species, use your own *.gmt file')
        } else if (input$useOwnPathways) {
            fileInput('gmtfile', 'Choose gmt file')
        } else if (detectedSpecies == 'mm') {
            includeHTML('static/mm-checkboxes.html')
        } else {
            includeHTML('static/hs-checkboxes.html')
        }
    })

    processButtonClick <- eventReactive(input$submitButton, {

        if (is.null(detectedSpecies))
            return(NULL)
        print(input$useOwnPathways)
        if (detectedSpecies == 'unknown') {
            # Using user defined pathways, therefore we don't convert input genes to entrez
            gmt.file <- input$gmtfile$datapath
            if (is.null(gmt.file))
                return(NULL)

            addClass(selector = "body", class = "sidebar-collapse")
            shinyjs::show("container")

            pathways <<- get_pathways(gmt.file)
            print(paste("User provided pathways:", length(pathways)))
        } else if (input$useOwnPathways) {
            # Using user defined pathways, therefore we don't convert input genes to entrez
            gmt.file <- input$gmtfile$datapath
            if (is.null(gmt.file))
                return(NULL)

            addClass(selector = "body", class = "sidebar-collapse")
            shinyjs::show("container")

            pathways <<- get_pathways(gmt.file)
            print(paste("User provided pathways:", length(pathways)))
        } else {
            print(paste(input$selectedGenesets))
            # Predefined pathways, converting genes to entrez
            if (length(input$selectedGenesets) == 0)
                return(NULL)

            addClass(selector = "body", class = "sidebar-collapse")
            shinyjs::show("container")

            convertToEntrez()
            pathways <<- NULL
            bonus_pathway_files <-
                sapply(input$selectedGenesets, function(set) { paste0('./data/', detectedSpecies, '.', set, '.gmt') })
            lapply(bonus_pathway_files, loadPathwayFile)

            print(paste("User selected predefined pathways:", length(pathways)))
        }

        res <- NULL
        set.seed(42)
        res <- fgsea(pathways, ranks, nperm=25000, maxSize=500)
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
