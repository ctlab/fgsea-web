
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
                  'Choose *.rnk file',
                  accept='.rnk'),
        fileInput('gmtfile',
                  'Choose *.gmt file',
                  accept='.gmt')
    ),
    dashboardBody(
        includeHTML("include.html"),
        dataTableOutput("contents")
    )
)
