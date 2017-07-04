
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinydashboard)
library(shinyjs)
library(DT)

SIDEBAR_WIDTH <- 320

dashboardPage(

    dashboardHeader(title = 'FGSEA',
                    titleWidth = SIDEBAR_WIDTH),

    dashboardSidebar(
        width = SIDEBAR_WIDTH,

        useShinyjs(),

        fileInput('tsvfile',
                  'Choose DE tsv file'),

        uiOutput("uiDEInfo"),
        uiOutput("selectPathways"),
        actionButton("submitButton", "Perform Enrichment", class="btn-submit"),

        downloadLink('downloadExampleTSV', 'Download sample DE tsv file', class="sample-download")
    ),

    dashboardBody(
        useShinyjs(),
        includeCSS("static/styles.css"),
        shinyjs::hidden(div(id = "container", div(id = "loader"))),
        includeHTML("static/include.html"),
        dataTableOutput("contents")
    )
)
