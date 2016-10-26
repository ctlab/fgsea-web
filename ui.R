
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyBS)
library(DT)

dashboardPage(

    dashboardHeader(title = 'FGSEA'),

    dashboardSidebar(

        useShinyjs(),
        extendShinyjs("app.js"),

        fileInput('rnkfile',
                  'Choose rnk file',
                  accept='text/tab-separated-values'),

        fileInput('gmtfile',
                  'Choose gmt file',
                  accept='text/tab-separated-values'),

        downloadButton('downloadExampleRNK', 'Download sample rnk file', class="sample-download"),
        downloadButton('downloadExampleGMT', 'Download sample gmt file', class="sample-download")
    ),

    dashboardBody(
        tags$head(
            tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
        ),
        includeHTML("include.html"),
        dataTableOutput("contents")
    )
)
