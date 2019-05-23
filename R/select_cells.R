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
    shiny::titlePanel("Choose cells for a subset"),

    # Sidebar layout with input and output definitions ----
    shiny::sidebarLayout(

      # Sidebar panel for inputs ----
      shiny::sidebarPanel(
        # done button
        shiny::actionButton("choose_toggle", "Choose/unchoose"),
        # clear button
        shiny::actionButton("reset", "Clear"),
        # done button
        shiny::actionButton("done", "Done"),
        h3("Instructions:"),
        tags$ol(
          tags$li("Highlight points by clicking and dragging."),
          tags$li("Click the 'Choose/unchoose' button."),
          tags$li("Repeat until all of the desired cells are black."),
          tags$li("Click 'Done'.")
        ),
        h4("Details:"),
        tags$ul(
          tags$li("To start over, click 'Clear'"),
          tags$li("You can also choose/unchoose specific cells by clicking on them directly")
        )
      ),

      # Main panel for displaying outputs ----
      shiny::mainPanel(
        shiny::plotOutput("plot1", height="auto",
                   click = "plot1_click",
                   brush = shiny::brushOpts(id = "plot1_brush"))
      )
    )
  )

  server <- function(input, output, session) {

    vals <- reactiveValues(
      keeprows = rep(FALSE, nrow(colData(cds)))
    )

    output$plot1 <- shiny::renderPlot({
      # Plot the kept and excluded points as two separate data sets
      colData(cds)$keep <- vals$keeprows

      suppressMessages(plot_cells(cds, cell_size = 1, label_cell_groups = FALSE) + geom_point(alpha = colData(cds)$keep))
    }, height = function() {
      session$clientData$output_plot1_width
    })

    # Toggle points that are clicked
    shiny::observeEvent(input$plot1_click, {
      res <- shiny::nearPoints(as.data.frame(reducedDims(cds)[[reduction_method]]), xvar = "V1", yvar = "V2", input$plot1_click, allRows = TRUE)
      vals$keeprows <- vals$keeprows | res$selected_
    })

    # Toggle points that are brushed, when button is clicked
    shiny::observeEvent(input$choose_toggle, {
      res <- shiny::brushedPoints(as.data.frame(reducedDims(cds)[[reduction_method]]), xvar = "V1", yvar = "V2", input$plot1_brush, allRows = TRUE)
      vals$keeprows <- vals$keeprows | res$selected_
    })

    # Reset all points
    shiny::observeEvent(input$reset, {
      vals$keeprows <- rep(FALSE, nrow(colData(cds)))
    })

    shiny::observeEvent(input$done, {
      shiny::stopApp(vals$keeprows)
    })

  }
  sel <- suppressMessages(shiny::runApp(shiny::shinyApp(ui, server)))
  if(return_list) {
    return(row.names(colData(cds)[sel,]))
  } else {
    return(cds[,sel])
  }
}


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
choose_graph_segment <- function(cds, reduction_method = "UMAP", return_list = FALSE) {
  ui <- shiny::fluidPage(
    shiny::titlePanel("Choose cells for a subset"),

    # Sidebar layout with input and output definitions ----
    shiny::sidebarLayout(

      # Sidebar panel for inputs ----
      shiny::sidebarPanel(
        # done button
        shiny::actionButton("choose_toggle", "Choose/unchoose"),
        # clear button
        shiny::actionButton("reset", "Clear"),
        # done button
        shiny::actionButton("done", "Done"),
        h3("Instructions:"),
        tags$ol(
          tags$li("Highlight points by clicking and dragging."),
          tags$li("Click the 'Choose/unchoose' button."),
          tags$li("Repeat until all of the desired cells are black."),
          tags$li("Click 'Done'.")
        ),
        h4("Details:"),
        tags$ul(
          tags$li("To start over, click 'Clear'"),
          tags$li("You can also choose/unchoose specific cells by clicking on them directly")
        )
      ),

      # Main panel for displaying outputs ----
      shiny::mainPanel(
        shiny::plotOutput("plot1", height="auto",
                          click = "plot1_click",
                          brush = shiny::brushOpts(id = "plot1_brush"))
      )
    )
  )

  server <- function(input, output, session) {

    vals <- reactiveValues(
      keeprows = rep(FALSE, nrow(colData(cds)))
    )

    output$plot1 <- shiny::renderPlot({
      # Plot the kept and excluded points as two separate data sets
      colData(cds)$keep <- vals$keeprows

      suppressMessages(plot_cells(cds, cell_size = 1, label_cell_groups = FALSE) + geom_point(alpha = colData(cds)$keep))
    }, height = function() {
      session$clientData$output_plot1_width
    })

    # Toggle points that are clicked
    shiny::observeEvent(input$plot1_click, {
      res <- shiny::nearPoints(as.data.frame(reducedDims(cds)[[reduction_method]]), xvar = "V1", yvar = "V2", input$plot1_click, allRows = TRUE)
      vals$keeprows <- vals$keeprows | res$selected_
    })

    # Toggle points that are brushed, when button is clicked
    shiny::observeEvent(input$choose_toggle, {
      res <- shiny::brushedPoints(as.data.frame(reducedDims(cds)[[reduction_method]]), xvar = "V1", yvar = "V2", input$plot1_brush, allRows = TRUE)
      vals$keeprows <- vals$keeprows | res$selected_
    })

    # Reset all points
    shiny::observeEvent(input$reset, {
      vals$keeprows <- rep(FALSE, nrow(colData(cds)))
    })

    shiny::observeEvent(input$done, {
      shiny::stopApp(vals$keeprows)
    })

  }
  sel <- suppressMessages(shiny::runApp(shiny::shinyApp(ui, server)))
  if(return_list) {
    return(row.names(colData(cds)[sel,]))
  } else {
    return(cds[,sel])
  }
}

