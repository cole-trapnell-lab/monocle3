
#' Orders cells according to pseudotime.
#'
#' Assigns cells a pseudotime value based on their projection on the principal
#' graph learned in the \code{learn_graph} function and the position of chosen
#' root states. This function takes as input a cell_data_set and returns it
#' with two new colData columns: \code{Pseudotime} and \code{State}, which
#' together encode where each cell maps to the principal graph trajectory.
#' \code{order_cells()} optionally takes "root" state(s) in the form of cell
#' or principal graph node IDs, which you can use to specify the start of the
#' trajectory. If you don't provide a root state, an plot will be generated
#' where you can choose the root state(s) interactively. The trajectory will be
#' composed of segments. The cells from a segment will share the same value of
#' \code{State}.
#'
#' @param cds the cell_data_set upon which to perform this operation
#' @param reduction_method a string specifying the reduced dimension method to
#'   use when ordering cells. Currently only "UMAP" is supported.
#' @param root_pr_nodes NULL or a vector of starting principal points. If
#'   provided, Pseudotime will start (i.e. be zero) at these graph nodes. Both
#'   \code{root_pr_nodes} and \code{root_cells} cannot be provided.
#' @param root_cells NULL or a vector of starting cells. If provided,
#'   Pseudotime will start (i.e. be zero) at these cells. Both
#'   \code{root_pr_nodes} and \code{root_cells} cannot be provided.
#' @param verbose Whether to show running information for order_cells
#'
#' @return an updated cell_data_set object, in which colData contains values
#'   for State and Pseudotime for each cell.
#' @export
order_cells <- function(cds,
                        reduction_method = "UMAP",
                        root_pr_nodes=NULL,
                        root_cells=NULL,
                        verbose = FALSE){

  assertthat::assert_that(is(cds, "cell_data_set"))
  assertthat::assert_that(assertthat::are_equal("UMAP", reduction_method),
                          msg = paste("Currently only 'UMAP' is accepted as a",
                                      "reduction_method."))
  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste0("No dimensionality reduction for ",
                                      reduction_method, " calculated. ",
                                      "Please run reduce_dimensions with ",
                                      "reduction_method = ", reduction_method,
                                      ", cluster_cells, and learn_graph ",
                                      "before running order_cells."))
  assertthat::assert_that(!is.null(cds@clusters[[reduction_method]]),
                          msg = paste("No cell clusters for",
                                      reduction_method, "calculated.",
                                      "Please run cluster_cells with",
                                      "reduction_method =", reduction_method,
                                      "and run learn_graph before running",
                                      "order_cells."))
  assertthat::assert_that(!is.null(principal_graph(cds)[[reduction_method]]),
                          msg = paste("No principal graph for",
                                      reduction_method, "calculated.",
                                      "Please run learn_graph with",
                                      "reduction_method =", reduction_method,
                                      "before running order_cells."))
  assertthat::assert_that(
    igraph::vcount(principal_graph(cds)[[reduction_method]]) < 10000,
    msg = paste("principal graph is too large. order_cells doesn't support",
                "more than 10 thousand centroids."))
  if(!is.null(root_pr_nodes)) {
    assertthat::assert_that(
      all(root_pr_nodes %in%
            igraph::V(principal_graph(cds)[[reduction_method]])$name),
      msg = paste("All provided root_pr_nodes must be present in the",
                  "principal graph."))
  }

  if(!is.null(root_cells)) {
    assertthat::assert_that(all(root_cells %in% row.names(colData(cds))),
                            msg = paste("All provided root_cells must be",
                                        "present in the cell data set."))
  }
  if(is.null(root_cells) & is.null(root_pr_nodes)) {
    assertthat::assert_that(interactive(),
                            msg = paste("When not in interactive mode, either",
                                        "root_pr_nodes or root_cells must be",
                                        "provided."))
  }
  assertthat::assert_that(!all(c(!is.null(root_cells),
                                 !is.null(root_pr_nodes))),
                            msg = paste("Please specify either root_pr_nodes",
                                        "or root_cells, not both."))

  if (is.null(root_pr_nodes) & is.null(root_cells)){
    if (interactive()){
      root_pr_nodes <-
        select_trajectory_roots(cds, reduction_method = reduction_method)
    }
  } else if(!is.null(root_cells)){
    closest_vertex <- cds@principal_graph_aux[[
      reduction_method]]$pr_graph_cell_proj_closest_vertex
    root_pr_nodes <- paste("Y_", closest_vertex[root_cells,], sep="")
  }

  cds@principal_graph_aux[[reduction_method]]$root_pr_nodes <- root_pr_nodes

  cc_ordering <- extract_general_graph_ordering(cds, root_pr_nodes, verbose,
                                                reduction_method)
  colData(cds)$Pseudotime <- cc_ordering[row.names(colData(cds)), ]$pseudo_time

  cds
}

extract_general_graph_ordering <- function(cds,
                                           root_cell,
                                           verbose=T,
                                           reduction_method) {
  Z <- t(reducedDims(cds)[[reduction_method]])
  Y <- cds@principal_graph_aux[[reduction_method]]$dp_mst
  pr_graph <- principal_graph(cds)[[reduction_method]]

  res <- list(subtree = pr_graph, root = root_cell)

  parents <- rep(NA, length(igraph::V(pr_graph)))
  states <- rep(NA, length(igraph::V(pr_graph)))

  if(any(is.na(igraph::E(pr_graph)$weight))) {
    igraph::E(pr_graph)$weight <- 1
  }

  # do pseudotime calculation on the cell-wise graph
  # 1. identify nearest cells to the selected principal node
  # 2. build a cell-wise graph for each louvain group
  # 3. run the distance function to assign pseudotime for each cell
  closest_vertex <- find_nearest_vertex(Y[, root_cell, drop = F], Z)
  closest_vertex_id <- colnames(cds)[closest_vertex]

  cell_wise_graph <-
    cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_tree
  cell_wise_distances <- igraph::distances(cell_wise_graph,
                                           v = closest_vertex_id)

  if (length(closest_vertex_id) > 1){
    node_names <- colnames(cell_wise_distances)
    pseudotimes <- apply(cell_wise_distances, 2, min)
  }else{
    node_names <- names(cell_wise_distances)
    pseudotimes <- cell_wise_distances
  }

  names(pseudotimes) <- node_names

  ordering_df <- data.frame(sample_name = igraph::V(cell_wise_graph)$name,
                            pseudo_time = as.vector(pseudotimes)
  )
  row.names(ordering_df) <- ordering_df$sample_name
  return(ordering_df)
}

# Select the roots of the principal graph
select_trajectory_roots <- function(cds, x=1, y=2, # nocov start
                                    reduction_method) {
  reduced_dim_coords <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst)

  ica_space_df <- as.data.frame(reduced_dim_coords)
  use_3d <- ncol(ica_space_df) >= 3
  if (use_3d){
    colnames(ica_space_df) = c("prin_graph_dim_1", "prin_graph_dim_2",
                               "prin_graph_dim_3")
  }
  else{
    colnames(ica_space_df) = c("prin_graph_dim_1", "prin_graph_dim_2")
  }

  ica_space_df$sample_name <- row.names(ica_space_df)
  ica_space_df$sample_state <- row.names(ica_space_df)

  dp_mst <- principal_graph(cds)[[reduction_method]]

  if (is.null(dp_mst)){
    stop("You must first call orderCells() before using this function")
  }

  if (use_3d){
    edge_df <- dp_mst %>%
      igraph::as_data_frame() %>%
      select_(source = "from", target = "to") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select_(source="sample_name",
                                    source_prin_graph_dim_1="prin_graph_dim_1",
                                    source_prin_graph_dim_2="prin_graph_dim_2",
                                    source_prin_graph_dim_3="prin_graph_dim_3"),
                       by = "source") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select_(target="sample_name",
                                    target_prin_graph_dim_1="prin_graph_dim_1",
                                    target_prin_graph_dim_2="prin_graph_dim_2",
                                    target_prin_graph_dim_3="prin_graph_dim_3"),
                       by = "target")
  }else{
    edge_df <- dp_mst %>%
      igraph::as_data_frame() %>%
      dplyr::select_(source = "from", target = "to") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select_(source="sample_name",
                                    source_prin_graph_dim_1="prin_graph_dim_1",
                                    source_prin_graph_dim_2="prin_graph_dim_2"),
                       by = "source") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select_(target="sample_name",
                                    target_prin_graph_dim_1="prin_graph_dim_1",
                                    target_prin_graph_dim_2="prin_graph_dim_2"),
                       by = "target")
  }

  num_roots <- nrow(ica_space_df)
  sel <- rep(FALSE, nrow(ica_space_df))

  if (use_3d){
    ui <- shiny::fluidPage(
      shiny::titlePanel("Choose your root nodes"),

      # Sidebar layout with input and output definitions ----
      shiny::sidebarLayout(

        # Sidebar panel for inputs ----
        shiny::sidebarPanel(
          # clear button
          shiny::actionButton("reset", "Clear"),
          # done button
          shiny::actionButton("done", "Done")
        ),

        # Main panel for displaying outputs ----
        shiny::mainPanel(
          plotlyOutput("plot1")
        )
      )
    )

    server <- function(input, output) {

      vals <- shiny::reactiveValues(
        keeprows = rep(TRUE, nrow(ica_space_df))
      )

      output$plot1 <- plotly::renderPlotly({
        ica_space_df$keep <- FALSE
        ica_space_df[ vals$keeprows,]$keep <- TRUE

        plotly::plot_ly(x = ica_space_df$prin_graph_dim_1,
                y = ica_space_df$prin_graph_dim_2,
                z = ica_space_df$prin_graph_dim_3,
                colors = c('blue', 'black'), key = ica_space_df$sample_name,
                color = ica_space_df$keep, size = 1,
                type = "scatter3d", mode="markers") %>%
          plotly::layout(showlegend = FALSE)
      })
      # Toggle points that are clicked
      shiny::observeEvent(event_data("plotly_click"), {
        d <- event_data("plotly_click")
        new_keep <- rep(FALSE, nrow(ica_space_df))
        new_keep[which(ica_space_df$sample_name == d$key)] <- TRUE
        vals$keeprows <- xor(vals$keeprows, new_keep)
      })

      # Reset all points
      shiny::observeEvent(input$reset, {
        vals$keeprows <- rep(TRUE, nrow(ica_space_df))
      })

      shiny::observeEvent(input$done, {
        stopApp(vals$keeprows)
      })

    }
    sel <- shiny::runApp(shiny::shinyApp(ui, server))
  } else {
    ui <- shiny::fluidPage(
      shiny::titlePanel("Choose your root nodes"),

      # Sidebar layout with input and output definitions ----
      shiny::sidebarLayout(

        # Sidebar panel for inputs ----
        shiny::sidebarPanel(
          # done button
          shiny::actionButton("choose_toggle", "Choose/unchoose"),
          # clear button
          shiny::actionButton("reset", "Clear"),
          # done button
          shiny::actionButton("done", "Done")
        ),

        # Main panel for displaying outputs ----
        shiny::mainPanel(
          shiny::plotOutput("plot1", height = 350,
                     click = "plot1_click",
                     brush = shiny::brushOpts(id = "plot1_brush"))
        )
      )
    )

    server <- function(input, output) {

      vals <- shiny::reactiveValues(
        keeprows = rep(TRUE, nrow(ica_space_df))
      )

      output$plot1 <- shiny::renderPlot({
        # Plot the kept and excluded points as two separate data sets
        keep    <- ica_space_df[ vals$keeprows, , drop = FALSE]
        exclude <- ica_space_df[!vals$keeprows, , drop = FALSE]

        ggplot(keep, aes(prin_graph_dim_1, prin_graph_dim_2)) +
          geom_point(alpha = .7) +
          geom_point(data = exclude, shape = 21, fill = NA, color = "blue") +
          geom_segment(data = edge_df,  aes(x = source_prin_graph_dim_1,
                                            xend = target_prin_graph_dim_1,
                                            y = source_prin_graph_dim_2,
                                            yend = target_prin_graph_dim_2)) +
          labs(x="Component 1", y="Component 2") +
          monocle3:::monocle_theme_opts()
      })

      # Toggle points that are clicked
      shiny::observeEvent(input$plot1_click, {
        res <- shiny::nearPoints(ica_space_df, input$plot1_click,
                                 allRows = TRUE)

        vals$keeprows <- xor(vals$keeprows, res$selected_)
      })

      # Toggle points that are brushed, when button is clicked
      shiny::observeEvent(input$choose_toggle, {
        res <- shiny::brushedPoints(ica_space_df, input$plot1_brush,
                                    allRows = TRUE)

        vals$keeprows <- xor(vals$keeprows, res$selected_)
      })

      # Reset all points
      shiny::observeEvent(input$reset, {
        vals$keeprows <- rep(TRUE, nrow(ica_space_df))
      })

      shiny::observeEvent(input$done, {
        stopApp(vals$keeprows)
      })

    }
    sel <- shiny::runApp(shiny::shinyApp(ui, server))
  }
  ## return indices of selected points
  as.character(ica_space_df$sample_name[which(!sel)])
} # nocov end

branch_nodes <- function(cds,reduction_method="UMAP"){
  g = principal_graph(cds)[[reduction_method]]
  branch_points <- which(igraph::degree(g) > 2)
  branch_points = branch_points[branch_points %in% root_nodes(cds) == FALSE]
  return(branch_points)
}

leaf_nodes <- function(cds,reduction_method="UMAP"){
  g = principal_graph(cds)[[reduction_method]]
  leaves <- which(igraph::degree(g) == 1)
  leaves = leaves[leaves %in% root_nodes(cds) == FALSE]
  return(leaves)
}

root_nodes <- function(cds, reduction_method="UMAP"){
  g = principal_graph(cds)[[reduction_method]]
  root_pr_nodes <- which(names(igraph::V(g)) %in%
                    cds@principal_graph_aux[[reduction_method]]$root_pr_nodes)
  names(root_pr_nodes) <-
    cds@principal_graph_aux[[reduction_method]]$root_pr_nodes
  return(root_pr_nodes)
}




