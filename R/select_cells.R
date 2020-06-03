#' Choose cells interactively to subset a cds
#'
#' @param cds CDS object to subset
#' @param reduction_method The reduction method to plot while choosing cells.
#' @param return_list Logical, return a list of cells instead of a subsetted
#'   CDS object.
#'
#' @return A subset CDS object. If return_list = FALSE, a list of cell names.
#' @export
#'
choose_cells <- function(cds,
                         reduction_method = c("UMAP", "tSNE", "PCA", "Aligned"),
                         return_list = FALSE) {
  reduction_method <- match.arg(reduction_method)
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste0("No dimensionality reduction for ",
                                       reduction_method, " calculated. ",
                                       "Please run reduce_dimensions with ",
                                       "reduction_method = ", reduction_method,
                                       ", cluster_cells, and learn_graph ",
                                       "before running choose_cells"))
  assertthat::assert_that(is.logical(return_list))
  assertthat::assert_that(interactive(),
                          msg = paste("choose_cells only works in",
                                      "interactive mode."))


  reduced_dims <- as.data.frame(reducedDims(cds)[[reduction_method]])
  names(reduced_dims)[1:2] <- c("V1", "V2")

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
        shiny::h3("Instructions:"),
        shiny::tags$ol(
          shiny::tags$li("Highlight points by clicking and dragging."),
          shiny::tags$li("Click the 'Choose/unchoose' button."),
          shiny::tags$li("Repeat until all of the desired cells are black."),
          shiny::tags$li("Click 'Done'.")
        ),
        shiny::h4("Details:"),
        shiny::tags$ul(
          shiny::tags$li("To start over, click 'Clear'"),
          shiny::tags$li(paste("You can also choose/unchoose specific cells",
                               "by clicking on them directly"))
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

    vals <- shiny::reactiveValues(
      keeprows = rep(FALSE, nrow(colData(cds)))
    )

    output$plot1 <- shiny::renderPlot({
      # Plot the kept and excluded points as two separate data sets
      colData(cds)$keep <- vals$keeprows

      suppressMessages(plot_cells(cds, reduction_method = reduction_method,
                                  cell_size = 1, label_cell_groups = FALSE,
                                  rasterize=FALSE) +
                         geom_point(alpha = colData(cds)$keep)) +
        theme(legend.position = "none")
    }, height = function() {
      session$clientData$output_plot1_width
    })

    # Toggle points that are clicked
    shiny::observeEvent(input$plot1_click, {
      res <- shiny::nearPoints(reduced_dims,
                               xvar = "V1", yvar = "V2", input$plot1_click,
                               allRows = TRUE)
      vals$keeprows <- vals$keeprows | res$selected_
    })

    # Toggle points that are brushed, when button is clicked
    shiny::observeEvent(input$choose_toggle, {
      res <- shiny::brushedPoints(reduced_dims,
                                  xvar = "V1", yvar = "V2", input$plot1_brush,
                                  allRows = TRUE)
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


#' Choose cells interactively along the path of a principal graph
#'
#' @param cds CDS object to be subsetted.
#' @param reduction_method The reduction method to plot while choosing cells.
#'   Currently only "UMAP" is supported.
#' @param return_list Logical, return a list of cells instead of a subsetted
#'   CDS object.
#' @param clear_cds Logical, clear CDS slots before returning.
#'   After clearing the cds, re-run processing from preprocess_cds(), ...
#'   Default is TRUE.
#'
#' @return A subset CDS object. If return_list = FALSE, a list of cell and
#'   graph node names.
#' @export
#'
choose_graph_segments <- function(cds,
                                 reduction_method = "UMAP",
                                 return_list = FALSE,
                                 clear_cds = TRUE) {

  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(assertthat::are_equal("UMAP", reduction_method),
                          msg = paste("Currently only 'UMAP' is accepted as a",
                                      "reduction_method."))
  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste0("No dimensionality reduction for ",
                                       reduction_method, " calculated. ",
                                       "Please run reduce_dimensions with ",
                                       "reduction_method = ", reduction_method,
                                       ", cluster_cells, and learn_graph ",
                                       "before running choose_graph_segments."))
  assertthat::assert_that(!is.null(cds@clusters[[reduction_method]]),
                          msg = paste("No cell clusters for",
                                      reduction_method, "calculated.",
                                      "Please run cluster_cells with",
                                      "reduction_method =", reduction_method,
                                      "and run learn_graph before running",
                                      "choose_graph_segments."))
  assertthat::assert_that(!is.null(principal_graph(cds)[[reduction_method]]),
                          msg = paste("No principal graph for",
                                      reduction_method, "calculated.",
                                      "Please run learn_graph with",
                                      "reduction_method =", reduction_method,
                                      "before running choose_graph_segments."))
  assertthat::assert_that(is.logical(return_list))
  assertthat::assert_that(interactive(),
                          msg = paste("choose_graph_segments only works in",
                                      "interactive mode."))

  dp_mst <- cds@principal_graph[[reduction_method]]

  princ_points <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>%
    as.data.frame() %>%
    dplyr::select_(x = 1, y = 2) %>%
    dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))
  row.names(princ_points) <- princ_points$sample_name

  edge_df <- dp_mst %>%
    igraph::as_data_frame() %>%
    dplyr::select_(source = "from", target = "to") %>%
    dplyr::left_join(princ_points %>%
                       dplyr::select_(source="sample_name",
                                      source_prin_graph_dim_1="x",
                                      source_prin_graph_dim_2="y"),
                     by = "source") %>%
    dplyr::left_join(princ_points %>%
                       dplyr::select_(target="sample_name",
                                      target_prin_graph_dim_1="x",
                                      target_prin_graph_dim_2="y"),
                     by = "target")

  data_df <- data.frame(reducedDims(cds)[[reduction_method]])

  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)

  data_df <- as.data.frame(cbind(data_df, colData(cds)))

  data_df$cell_color = tryCatch({
    partitions(cds, reduction_method = reduction_method)[data_df$sample_name]},
    error = function(e) {NULL})
  data_df$chosen_cells <- FALSE

  ui <- shiny::fluidPage(
    shiny::titlePanel("Choose cells along a graph path"),

    # Sidebar layout with input and output definitions ----
    shiny::sidebarLayout(

      # Sidebar panel for inputs ----
      shiny::sidebarPanel(
        shiny::actionButton("choose_start", "Choose starting node"),
        shiny::actionButton("choose_end", "Choose ending nodes"),
        shiny::actionButton("connect", "Connect nodes"),
        # clear button
        shiny::actionButton("reset", "Clear"),
        # done button
        shiny::actionButton("done", "Done"),
        shiny::h3("Instructions:"),
        shiny::tags$ol(
          shiny::tags$li("Highlight starting principal_graph node."),
          shiny::tags$li("Click 'Choose starting node' to highlight."),
          shiny::tags$li("Highlight ending principal_graph node(s)."),
          shiny::tags$li("Click 'Choose ending nodes' to highlight."),
          shiny::tags$li(paste("Click 'Connect nodes' to highlight connecting nodes",
                        "and cells")),
          shiny::tags$li("Click 'Done' to return the chosen subset.")
        ),
        shiny::h4("Details:"),
        shiny::tags$ul(
          shiny::tags$li("To start over, click 'Clear'"),
          shiny::tags$li(paste("You can choose multiple ending nodes, but only 1",
                        "starting node."))
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

    vals <- shiny::reactiveValues(
      start = rep(FALSE, nrow(princ_points)),
      end =  rep(FALSE, nrow(princ_points)),
      chosen = rep(FALSE, nrow(princ_points)),
      chosen_cells = rep(FALSE, nrow(data_df))
    )

    output$plot1 <- shiny::renderPlot({
      # Plot the kept and excluded points as two separate data sets
      princ_points$start <- vals$start
      princ_points$end <- vals$end
      princ_points$chosen <- "Unchosen"
      princ_points$chosen[vals$start] <- "Start"
      princ_points$chosen[vals$end] <- "End"
      princ_points$chosen[vals$chosen] <- "Chosen"
      data_df$chosen_cells <- "gray"
      data_df$chosen_cells[vals$chosen_cells] <- "purple"
      suppressMessages(plot_principal_graph(cds, data_df, princ_points,
                                            label_branch_points = FALSE,
                                            label_leaves = FALSE,
                                            label_roots = FALSE))
    }, height = function() {
      session$clientData$output_plot1_width
    })

    # Toggle points that are brushed, when button is clicked
    shiny::observeEvent(input$choose_start, {
      res <- shiny::brushedPoints(princ_points, xvar = "x", yvar = "y",
                                  input$plot1_brush, allRows = TRUE)
      vals$start <- res$selected_
    })

    # Toggle points that are brushed, when button is clicked
    shiny::observeEvent(input$choose_end, {
      res <- shiny::brushedPoints(princ_points, xvar = "x", yvar = "y",
                                  input$plot1_brush, allRows = TRUE)
      vals$end <- vals$end | res$selected_
    })

    shiny::observeEvent(input$connect, {
      chosen <- tryCatch(
        get_principal_path(cds, reduction_method,
                           starting_cell = row.names(princ_points)[vals$start],
                           end_cells = row.names(princ_points)[vals$end]),
        error = function(e) print(e))
      vals$chosen <- vals$chosen | row.names(princ_points) %in% chosen$nodes
      vals$chosen_cells <- vals$chosen_cells | row.names(pData(cds)) %in%
        chosen$cells
      vals$start = rep(FALSE, nrow(princ_points))
      vals$end =  rep(FALSE, nrow(princ_points))
    })

    # Reset all points
    shiny::observeEvent(input$reset, {
      vals$start = rep(FALSE, nrow(princ_points))
      vals$end =  rep(FALSE, nrow(princ_points))
      vals$chosen =  rep(FALSE, nrow(princ_points))
      vals$chosen_cells = rep(FALSE, nrow(data_df))
    })

    shiny::observeEvent(input$done, {
      shiny::stopApp(list(nodes = row.names(princ_points)[vals$chosen],
                          cells = row.names(data_df)[vals$chosen_cells]))
    })

  }

  sel <- suppressMessages(shiny::runApp(shiny::shinyApp(ui, server)))

  if(return_list) {
    return(sel)
  } else {
    cds<-cds[,sel$cells]
    if( clear_cds )
      cds<-clear_cds_slots(cds)
    return(cds)
  }
}


get_principal_path <- function(cds, reduction_method,
                               starting_cell, end_cells) {

  subset_principal_nodes <- c()
  dp_mst <- principal_graph(cds)[[reduction_method]]

  for(end_cell in end_cells) {
    traverse_res <- traverse_graph(dp_mst, starting_cell, end_cell)
    path_cells <- names(traverse_res$shortest_path[[1]])
    if(length(path_cells) == 0) {
      stop(paste0("Starting and ending nodes are not connected"))
    }

    subset_principal_nodes <- c(subset_principal_nodes, path_cells)
  }
  subset_principal_nodes <- unique(subset_principal_nodes)
  corresponding_cells <- which(paste0("Y_", cds@principal_graph_aux[[
    reduction_method]]$pr_graph_cell_proj_closest_vertex) %in%
      subset_principal_nodes)
  subset_cells <- row.names(cds@principal_graph_aux[[
    reduction_method]]$pr_graph_cell_proj_closest_vertex)[corresponding_cells]

  return(list(nodes = subset_principal_nodes, cells =subset_cells))
}

traverse_graph <- function(g, starting_cell, end_cells){
  distance <- igraph::shortest.paths(g, v=starting_cell, to=end_cells)
  branchPoints <- which(igraph::degree(g) == 3)
  path <- igraph::shortest_paths(g, from = starting_cell, end_cells)

  return(list(shortest_path = path$vpath, distance = distance,
              branch_points = intersect(branchPoints, unlist(path$vpath))))
}


plot_principal_graph <- function(cds,
                                 data_df,
                                 princ_points,
                                 reduction_method = "UMAP",
                                 trajectory_graph_color="black",
                                 trajectory_graph_segment_size=0.75,
                                 label_groups_by_cluster=TRUE,
                                 group_label_size=2,
                                 labels_per_group=1,
                                 label_branch_points=TRUE,
                                 label_roots=TRUE,
                                 label_leaves=TRUE,
                                 graph_label_size=2,
                                 cell_size=0.35,
                                 alpha = 1,
                                 min_expr=0.1,
                                 rasterize=FALSE) {

  gene_short_name <- NA
  sample_name <- NA
  #sample_state <- colData(cds)$State
  data_dim_1 <- NA
  data_dim_2 <- NA
  plotting_func <- ggplot2::geom_point


  ## Graph info

  ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>%
    as.data.frame() %>%
    dplyr::select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
    dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))

  dp_mst <- cds@principal_graph[[reduction_method]]

  edge_df <- dp_mst %>%
    igraph::as_data_frame() %>%
    dplyr::select_(source = "from", target = "to") %>%
    dplyr::left_join(ica_space_df %>%
                       dplyr::select_(
                         source="sample_name",
                         source_prin_graph_dim_1="prin_graph_dim_1",
                         source_prin_graph_dim_2="prin_graph_dim_2"),
                     by = "source") %>%
    dplyr::left_join(ica_space_df %>%
                       dplyr::select_(
                         target="sample_name",
                         target_prin_graph_dim_1="prin_graph_dim_1",
                         target_prin_graph_dim_2="prin_graph_dim_2"),
                     by = "target")

  g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2))

  g <- g + geom_point(color=data_df$chosen_cells, size=I(cell_size),
                      na.rm = TRUE, alpha = I(alpha))

  message(paste("cluster_cells() has not been called yet, can't color cells",
                "by cluster"))

  g <- g + geom_point(aes(x = x, y = y, color = chosen), data=princ_points) +
    scale_color_manual(values = c("Start" = "green", "End" = "blue",
                                  "Unchosen" = "black", "Chosen" = "purple"))
  g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1",
                                   y="source_prin_graph_dim_2",
                                   xend="target_prin_graph_dim_1",
                                   yend="target_prin_graph_dim_2"),
                        size=trajectory_graph_segment_size,
                        linetype="solid",
                        na.rm=TRUE,
                        data=edge_df)


  if (label_branch_points){
    mst_branch_nodes <- branch_nodes(cds)
    branch_point_df <- ica_space_df %>%
      dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
      dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))

    g <- g +
      geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                 shape = 21, stroke=I(trajectory_graph_segment_size),
                 color="white",
                 fill="black",
                 size=I(graph_label_size * 1.5),
                 na.rm=TRUE, branch_point_df) +

      geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                           label="branch_point_idx"),
                size=I(graph_label_size), color="white", na.rm=TRUE,
                branch_point_df)
  }

  if (label_leaves){
    mst_leaf_nodes <- leaf_nodes(cds)
    leaf_df <- ica_space_df %>%
      dplyr::slice(match(names(mst_leaf_nodes), sample_name)) %>%
      dplyr::mutate(leaf_idx = seq_len(dplyr::n()))

    g <- g +
      geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                 shape = 21, stroke=I(trajectory_graph_segment_size),
                 color="black",
                 fill="lightgray",
                 size=I(graph_label_size * 1.5),
                 na.rm=TRUE,
                 leaf_df) +
      geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                           label="leaf_idx"),
                size=I(graph_label_size), color="black", na.rm=TRUE, leaf_df)
  }

  if (label_roots){
    mst_root_nodes <- root_nodes(cds)
    root_df <- ica_space_df %>%
      dplyr::slice(match(names(mst_root_nodes), sample_name)) %>%
      dplyr::mutate(root_idx = seq_len(dplyr::n()))

    g <- g +
      geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                 shape = 21, stroke=I(trajectory_graph_segment_size),
                 color="black",
                 fill="white",
                 size=I(graph_label_size * 1.5),
                 na.rm=TRUE,
                 root_df) +
      geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                           label="root_idx"),
                size=I(graph_label_size), color="black", na.rm=TRUE, root_df)
  }

  g <- g +
    monocle_theme_opts() +
    xlab(paste(reduction_method, 1)) +
    ylab(paste(reduction_method, 2)) +
    theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill='white'))
  g
}



