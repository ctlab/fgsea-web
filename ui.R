
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinydashboard)
library(shinyjs)
library(DT)

dashboardPage(

    dashboardHeader(title = 'FGSEA'),

    dashboardSidebar(

        useShinyjs(),

        fileInput('rnkfile',
                  'Choose rnk file'),

        fileInput('gmtfile',
                  'Choose gmt file'),

        checkboxGroupInput("checkGroup", label = h3("Select gene sets"),
                           choices = list(
                               "Hallmark gene sets" = 'hs.hallmark',
                               "Positional gene sets" = 'hs.c1',
                               "Chemical and genetic perturbations" = 'hs.cgp',
                               #"Canonical pathways" = 'hs.cp',
                               "BioCarta gene sets" = 'hs.biocarta',
                               "KEGG gene sets" = 'hs.kegg',
                               "Reactome gene sets" = 'hs.reactome',
                               #"Motif gene sets" = 'hs.c3',
                               "MicroRNA targets" = 'hs.mir',
                               "Transcription factor targets" = 'hs.tft',
                               "Cancer gene neighborhoods" = 'hs.cgn',
                               "Cancer modules" = 'hs.cm',
                               "GO gene sets" = 'hs.c5',
                               "Oncogenic signatures" = 'hs.c6',
                               "immunologic signatures" = 'hs.c7'
                               )
                           ),

        downloadLink('downloadExampleRNK', 'Download sample rnk file', class="sample-download"),
        downloadLink('downloadExampleGMT', 'Download sample gmt file', class="sample-download")
    ),

    dashboardBody(
        useShinyjs(),
        tags$head(
            tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
        ),
        shinyjs::hidden(div(id = "container", div(id = "loader"))),
        includeHTML("include.html"),
        dataTableOutput("contents")
    )
)
