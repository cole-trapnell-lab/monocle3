#' Choose cells interactively to make subset cds
#'
#' @param cds CDS object
#' @param reduction_method
#' @param return_list Logical, return a list of cells.
#'
#' @return A subset CDS object. If return_list = FALSE, a list of cell names.
#' @export
#'
#' @examples
choose_cells <- function(cds, reduction_method = "UMAP", return_list = FALSE) {
  ui <- shiny::fluidPage(
    titlePanel("Choose cells for a subset"),

    # Sidebar layout with input and output definitions ----
    sidebarLayout(

      # Sidebar panel for inputs ----
      sidebarPanel(
        # done button
        actionButton("choose_toggle", "Choose/unchoose"),
        # clear button
        actionButton("reset", "Clear"),
        # done button
        actionButton("done", "Done")
      ),

      # Main panel for displaying outputs ----
      mainPanel(
        plotOutput("plot1", height="auto",
                   click = "plot1_click",
                   brush = brushOpts(id = "plot1_brush"))
      )
    )
  )

  server <- function(input, output, session) {

    vals <- reactiveValues(
      keeprows = rep(FALSE, nrow(colData(cds)))
    )

    output$plot1 <- renderPlot({
      # Plot the kept and excluded points as two separate data sets
      colData(cds)$keep <- vals$keeprows

      plot_cells(cds) + geom_point(alpha = colData(cds)$keep)
    }, height = function() {
      session$clientData$output_plot1_width
    })

    # Toggle points that are clicked
    observeEvent(input$plot1_click, {
      res <- nearPoints(as.data.frame(reducedDims(cds)[[reduction_method]]), xvar = "V1", yvar = "V2", input$plot1_click, allRows = TRUE)

      vals$keeprows <- xor(vals$keeprows, res$selected_)
    })

    # Toggle points that are brushed, when button is clicked
    observeEvent(input$choose_toggle, {
      res <- brushedPoints(as.data.frame(reducedDims(cds)[[reduction_method]]), xvar = "V1", yvar = "V2", input$plot1_brush, allRows = TRUE)

      vals$keeprows <- xor(vals$keeprows, res$selected_)
    })

    # Reset all points
    observeEvent(input$reset, {
      vals$keeprows <- rep(FALSE, nrow(colData(cds)))
    })

    observeEvent(input$done, {
      stopApp(vals$keeprows)
    })

  }
  sel <- runApp(shinyApp(ui, server))
  if(return_list) {
    return(row.names(colData(cds)[sel,]))
  } else {
    return(cds[,sel])
  }
}
