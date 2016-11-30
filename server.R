
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

pathways <- NULL
ranks <- NULL
stripped_to_full <- NULL
full_to_stripped <- NULL
species <- NULL

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

convertToEntrez <- function(genes, format) {
    require(org.Hs.eg.db)
    require(org.Mm.eg.db)

    converted <- AnnotationDbi::mapIds(org.Mm.eg.db, keys=genes, column="ENTREZID", keytype=format, multiVals="first")
    species <<- 'mm'
    if (!exists('converted')) {
        converted <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=genes, column="ENTREZID", keytype=format, multiVals="first")
        species <<- 'hs'
    }
    print(paste('failed to convert', sum(is.na(converted)), 'genes'))
    inds <- which(!is.na(converted))
    ranks <<- ranks[inds]
    names(ranks) <<- converted[inds]
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

    output$genesets <- renderUI({
        checkbox = NULL
        if (input$radio == 'mm') {
            checkboxGroupInput("selectedGenesets", label = h3("Select gene sets"),
                               choices = list(
                                   "H hallmark gene sets" = 'mm.hallmark',
                                   "C2 curated gene sets" = 'mm.c2',
                                   "C3 motif gene sets" = 'mm.c3',
                                   "C4 computational gene sets" = 'mm.c4',
                                   "C5 GO gene sets" = 'mm.c5',
                                   "C6 oncogenic signatures" = 'mm.c6',
                                   "C7 immunologic signatures" = 'mm.c7'
                               )
            )
        } else {
            checkboxGroupInput("selectedGenesets", label = h3("Select gene sets"),
                               choices = list(
                                   "H hallmark gene sets" = 'hs.hallmark',
                                   "C1 positional gene sets" = 'hs.c1',
                                   "C2 curated gene sets" = 'hs.c2',
                                   "C3 motif gene sets" = 'hs.c3',
                                   "C4 computational gene sets" = 'hs.c4',
                                   "C5 GO gene sets" = 'hs.c5',
                                   "C6 oncogenic signatures" = 'hs.c6',
                                   "C7 immunologic signatures" = 'hs.c7'
                               )
            )
        }

    })

    output$contents <- renderDataTable(
        {
            rnk.file <- input$rnkfile$datapath
            gmt.file <- input$gmtfile$datapath

            if (is.null(rnk.file) | is.null(gmt.file))
                return(NULL)

            # hide sidebar and show loader animation
            addClass(selector = "body", class = "sidebar-collapse")
            shinyjs::show("container")

            # load rnk file
            ranks <<- get_ranks(rnk.file)

            # convert to entrez
            idExample <- names(ranks)[1]
            if (is.na(as.numeric(idExample))) {
                if (startsWith(idExample, 'ENS')) {
                    print('detected format: ensembl')
                    convertToEntrez(names(ranks), "ENSEMBL")
                } else {
                    print('detected format: symbol')
                    convertToEntrez(names(ranks), "SYMBOL")
                }
            }


            # load pathways from user
            pathways <<- get_pathways(gmt.file)
            print(paste("User provided pathways:", length(pathways)))

            # load base pathways
            bonus_pathway_files <- sapply(input$selectedGenesets, function(set) { paste0('./data/', set, '.gmt') })
            if (length(bonus_pathway_files > 0)) {
                lapply(bonus_pathway_files, loadPathwayFile)
            }
            print(paste("Total pathways:", length(pathways)))

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
            order = list(list(2, 'asc'))
        )
    )

    observe({
        if (!is.null(input$pathway)) {
            pathway <- stripped_to_full[[input$pathway]]
            filePath <- tempfile(tmpdir = "www", fileext = 'plot.svg')
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
