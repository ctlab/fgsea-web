
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

        uiOutput("useOwnPathwaysRadio"),
        uiOutput("selectPathways"),
        actionButton("submitButton", "Perform Enrichment"),

        downloadLink('downloadExampleRNK', 'Download sample rnk file', class="sample-download"),
        downloadLink('downloadExampleGMT', 'Download sample gmt file', class="sample-download")
    ),

    dashboardBody(
        useShinyjs(),
        includeCSS("styles.css"),
        shinyjs::hidden(div(id = "container", div(id = "loader"))),
        includeHTML("include.html"),
        dataTableOutput("contents")
    )
)
