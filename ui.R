
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("Fast Gene Set Enrichment Analysis"),

  # Sidebar with a slider input for number of bins
  sidebarPanel(
      fileInput('rnkfile', 
                'Choose *.rnk file', 
                accept='.rnk'),
      fileInput('gmtfile', 
                'Choose *.gmt file', 
                accept='.gmt')
  ),

  mainPanel(
    tableOutput("contents")
  )
))
