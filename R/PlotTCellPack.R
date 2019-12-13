#' Create T cell pack
#'
#' Creates a T cell pack from a gliph file, single-cell gene expression, or a combination of gliph and single-cell gene expression or gliph and T cell receptor sequencing.
#'
#' @param gliph Character vector providing the path to GLIPH_TCR_Table-convergence-groups.txt.
#' @param clonotype.data Optional 2 column data frame named "clonotype" and "frequency" or "clonotype" and "data".  If the second column
#' is "frequency", then the T cell clone size is drawn proportional to its frequency.  If the second column is "data", then the T cell clone
#' is colored according the variable (only discrete variables are supported at this time).
#' @param cell.data Optional 3 column data frame named "clonotype", "cell", and "data".
#' @param specificity.color Character vector providing the color (single value) of the specificity circles.
#' @param clonotype.color Character vector providing the color (single value) of the clonotype circles.
#' @param cell.color Character vector providing the color (single value)  of the cell circles.
#' @param color.gradient Character vector providing the RColorBrewer color palette used color circles with discrete or continuous values.
#' @param line.color Character vector providing the color (single value) of the circle borders.
#' @param label.color Character vector providing the color (single value) of the labels.
#' @param label.size Integer corresponding to text size of the labels.
#' @param label Character vector indicating if labels should be displayed.  Options include "none", "specificity", "clonotype", "cell", and "data".
#' @param legend Logic indicating if legend should be displayed (TRUE) or not (FALSE).
#' @param cluster.size Minimum number of clonotypes per specificity group.
#'
#' @return Returns a T cell circle packing plot.
#'
#' @examples
#' # T cell clone size proportional to frequency with GLIPH specificity groups
#' PlotTCellPack(gliph = gliph.example, clonotype.data = clonotype.data.example, legend = TRUE)
#'
#' # T cell colored by continuous variable with GLIPH specificity groups
#' PlotTCellPack(gliph = gliph.example, cell.data = cell.data.continuous.example, legend = TRUE)
#'
#' # T cell colored by discrete variable with GLIPH specificity groups
#' PlotTCellPack(gliph = gliph.example, cell.data = cell.data.discrete.example, legend = TRUE)
#'
#' # T cell colored by discrete variable without GLIPH specificity groups
#' PlotTCellPack(cell.data = cell.data.discrete.example, legend = TRUE)
#'
#' # T cell colored by continuous variable without GLIPH specificity groups
#' PlotTCellPack(cell.data = cell.data.continuous.example, legend = TRUE)
#'
#' @export
#' @import ggplot2
#' @import ggraph
#' @import igraph
#' @import RColorBrewer
#' @importFrom plyr ldply
#' @importFrom stringr str_split
#' @importFrom data.table fread
#' @importFrom grDevices colorRampPalette

PlotTCellPack <- function(gliph = NULL,
                          clonotype.data = NULL,
                          cell.data = NULL,
                          specificity.color = "#08306b",
                          clonotype.color = "#4292c6",
                          cell.color = "#9ecae1",
                          color.gradient = "RdYlBu",
                          line.color = "black",
                          label.color = "black",
                          label.size = 1,
                          label = "none",
                          legend = FALSE,
                          cluster.size = 1){

  # Input validation
  if(!is.null(cell.data)){
    if(!(all(names(cell.data) %in% c("clonotype", "cell", "data")))) {stop("cell.data must have the columns 'clonotype', 'cell', and 'data'", call. = TRUE)}
    if(any(is.na(cell.data$data))) {warning("cell.data cannot contain NA", call. = FALSE)}
  }

  if(!is.null(clonotype.data)){
    if(!(all(names(clonotype.data) %in% c("clonotype", "frequency")) | all(names(clonotype.data) %in% c("clonotype", "frequency", "data")))) {warning("clonotype.data must have the columns 'clonotype' and 'frequency' or 'clonotype' and 'data'", call. = TRUE)}
    if(any(clonotype.data$frequency %in% 0 & is.na(clonotype.data$frequency))) {warning("clonotype.data cannot contain NA or 0 values", call. = TRUE)}
  }

  if(!is.null(gliph) & any(class(gliph) %in% "character")){
    gliph <- data.frame(data.table::fread(gliph))
  }

  if(!is.null(gliph) & any(class(gliph) %in% "data.table")){
    gliph <- data.frame(gliph)
  }

  # Reformat GLIPH data
  if(!is.null(gliph)) {
    specifities <- stringr::str_split(string = gliph[gliph[,1] > cluster.size, 3], pattern = " ")
    names(specifities) <- 1:length(specifities)
    specifities <- plyr::ldply(specifities, data.frame)
    names(specifities) <- c("specificity", "clonotype")
  }

  # Only GLIPH data provided
  if(!is.null(gliph) & is.null(clonotype.data) & is.null(cell.data)){
    edges <- data.frame(from = specifities$specificity, to = specifities$clonotype)

    # Graph object
    graph <- igraph::graph_from_data_frame(d = edges)

    # Plot
    plot <- ggraph::ggraph(graph, layout = "circlepack") +
      ggraph::geom_node_circle(ggplot2::aes(fill = factor(depth)), size = 0.25, n = 50, color = line.color) +
      ggplot2::scale_fill_manual(values = c("0" = specificity.color, "1" = clonotype.color),
                                 name = "", labels = c("Specificity", "Clonotype")) +
      ggplot2::theme_void() +
      ggplot2::coord_fixed()

    # Define depth level for optional label
    clonotype = 1
    specificity = 0
  }

  # Only GLIPH and clonotype frequency data provided
  if(!is.null(gliph) & !is.null(clonotype.data) & is.null(cell.data) & any(names(clonotype.data) %in% "frequency")){
    specifities <- merge(specifities, clonotype.data, all = FALSE)
    specifities.aggregate <- aggregate(data = specifities, frequency~specificity, sum)
    names(specifities.aggregate) <- c("name", "frequency")
    names(clonotype.data) <- c("name", "frequency")
    clonotype.data <- clonotype.data[clonotype.data$name %in% specifities$clonotype,]
    vertices <- rbind(clonotype.data, specifities.aggregate)
    edges <- data.frame(from = specifities$specificity, to = specifities$clonotype)

    # Graph object
    graph <- igraph::graph_from_data_frame(d = edges, vertices = vertices)

    # Plot
    plot <- ggraph::ggraph(graph, layout = "circlepack", weight = frequency) +
      ggraph::geom_node_circle(ggplot2::aes(fill = factor(depth)), size = 0.25, n = 50, color = line.color) +
      ggplot2::scale_fill_manual(values = c("0" = specificity.color, "1" = clonotype.color),
                                 name = "", labels = c("Specificity", "Clonotype")) +
      ggplot2::theme_void() +
      ggplot2::coord_fixed()

    # Define depth level for optional label
    clonotype = 1
    specificity = 0
  }

  # Only GLIPH and discrete clonotype data provided
  if(!is.null(gliph) & !is.null(clonotype.data) & is.null(cell.data) & any(names(clonotype.data) %in% "data")){
    specifities <- merge(specifities, clonotype.data, all = FALSE)
    edges <- data.frame(from = specifities$specificity, to = specifities$clonotype)
    vertices <- data.frame(name = specifities$clonotype, data = specifities$data)
    vertices <- rbind(vertices, unique(data.frame(name = specifities$specificity, data = rep("1", nrow(specifities)))))

    # Graph object
    graph <- igraph::graph_from_data_frame(d = edges, vertices = vertices)

    # Plot
    getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, color.gradient))
    plot <- ggraph::ggraph(graph, layout = "circlepack") +
      ggraph::geom_node_circle(ggplot2::aes(fill = data), size = 0.25, n = 50, color = line.color) +
      ggraph::geom_node_circle(ggplot2::aes(fill = data, filter = leaf), size = 0.25, n = 50, color = line.color) +
      ggplot2::scale_fill_manual(values = c(specificity.color, getPalette(length(unique(specifities$data)))), breaks = sort(unique(specifities$data))) +
      ggplot2::theme_void() +
      ggplot2::coord_fixed() +
      ggplot2::labs(fill = "")

    # Define depth level for optional label
    cell = NULL
    clonotype = 1
    specificity = 0
  }

  # Only GLIPH and continuous cell data provided
  if(!is.null(gliph) & !is.null(cell.data) & class(cell.data$data) %in% c("integer", "numeric")){
    specifities <- merge(specifities, cell.data, all = FALSE)
    edges <- data.frame(from = specifities$specificity, to = specifities$clonotype)
    edges <- rbind(edges, data.frame(from = specifities$clonotype, to = specifities$cell))
    specifities.aggregate <- aggregate(data = specifities, data~specificity, mean)
    names(specifities.aggregate) <- c("name", "data")
    clonotype.aggregate <- aggregate(data = specifities, data~clonotype, mean)
    names(clonotype.aggregate) <- c("name", "data")
    vertices <- do.call("rbind", list(data.frame(name = specifities$cell, data = specifities$data),
                                      specifities.aggregate, clonotype.aggregate))
    # Graph object
    graph <- igraph::graph_from_data_frame(d = edges, vertices = vertices)

    # Plot
    plot <- ggraph::ggraph(graph, layout = "circlepack") +
      ggraph::geom_node_circle(ggplot2::aes(fill = data), size = 0.25, n = 50, color = line.color) +
      ggplot2::scale_fill_distiller(palette = color.gradient) +
      ggplot2::theme_void() +
      ggplot2::coord_fixed() +
      ggplot2::labs(fill = "")

    # Define depth level for optional label
    cell = 2
    clonotype = 1
    specificity = 0
  }

  # Only GLIPH and discrete cell data provided
  if(!is.null(gliph) & !is.null(cell.data) & class(cell.data$data) %in% c("factor", "character")){
    specifities <- merge(specifities, cell.data, all = FALSE)
    edges <- data.frame(from = specifities$specificity, to = specifities$clonotype)
    edges <- rbind(edges, data.frame(from =  specifities$clonotype, to = specifities$cell))
    vertices <- data.frame(name = specifities$cell, data = specifities$data)
    vertices <- rbind(vertices, unique(data.frame(name = specifities$specificity, data = rep("1", nrow(specifities)))))
    vertices <- rbind(vertices, unique(data.frame(name = specifities$clonotype, data = rep("2", nrow(specifities)))))

    # Graph object
    graph <- igraph::graph_from_data_frame(d = edges, vertices = vertices)

    # Plot
    getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, color.gradient))
    plot <- ggraph::ggraph(graph, layout = "circlepack") +
      ggraph::geom_node_circle(ggplot2::aes(fill = data), size = 0.25, n = 50, color = line.color) +
      ggraph::geom_node_circle(ggplot2::aes(fill = data, filter = leaf), size = 0.25, n = 50, color = line.color) +
      ggplot2::scale_fill_manual(values = c(specificity.color, clonotype.color, getPalette(length(unique(specifities$data)))), breaks = sort(unique(specifities$data))) +
      ggplot2::theme_void() +
      ggplot2::coord_fixed() +
      ggplot2::labs(fill = "")

    # Define depth level for optional label
    cell = 2
    clonotype = 1
    specificity = 0
  }

# Only discrete cell data provided
 if(is.null(gliph) & !is.null(cell.data) & class(cell.data$data) %in% c("factor", "character")){
   edges <- data.frame(from = cell.data$clonotype, to = cell.data$cell)
   vertices <- unique(data.frame(name = cell.data$cell, data = cell.data$data))
   vertices <- rbind(vertices, unique(data.frame(name = cell.data$clonotype, data = rep("1", nrow(cell.data)))))

   # Graph object
   graph <- igraph::graph_from_data_frame(d = edges, vertices = vertices)

   # Plot
   getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, color.gradient))
   plot <- ggraph::ggraph(graph, layout = "circlepack") +
     ggraph::geom_node_circle(ggplot2::aes(fill = data), size = 0.25, n = 50, color = line.color) +
     ggraph::geom_node_circle(ggplot2::aes(fill = data, filter = leaf), size = 0.25, n = 50, color = line.color) +
     ggplot2::scale_fill_manual(values = c(clonotype.color, getPalette(length(unique(cell.data$data)))), breaks = sort(unique(cell.data$data))) +
     ggplot2::theme_void() +
     ggplot2::coord_fixed() +
     ggplot2::labs(fill = "")

   # Define depth level for optional label
   cell = 1
   clonotype = 0
 }

# Only continous cell data provided
if(is.null(gliph) & !is.null(cell.data) & class(cell.data$data) %in% c("integer", "numeric")){
  edges <- data.frame(from = cell.data$clonotype, to = cell.data$cell)

  vertices <- unique(data.frame(name = cell.data$cell, data = cell.data$data))
  clonotype.aggregate <- aggregate(data = cell.data, data~clonotype, mean)
  names(clonotype.aggregate) <- c("name", "data")
  vertices <- rbind(vertices, clonotype.aggregate)

  # Graph object
  graph <- igraph::graph_from_data_frame(d = edges, vertices = vertices)

  # Plot
  getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, color.gradient))
  plot <- ggraph::ggraph(graph, layout = "circlepack") +
    ggraph::geom_node_circle(ggplot2::aes(fill = data), size = 0.25, n = 50, color = line.color) +
    ggplot2::scale_fill_distiller(palette = color.gradient) +
    ggplot2::theme_void() +
    ggplot2::coord_fixed() +
    ggplot2::labs(fill = "")

  # Define depth level for optional label
  cell = 1
  clonotype = 0
}

  if(legend == FALSE){
    plot <- plot + ggplot2::theme(legend.position="none")
  }

  if(label == "cell"){
    plot <- plot + ggraph::geom_node_text(aes(label=name, filter = depth == cell), color = label.color, size = label.size)
  }

  if(label == "data"){
    if(is.null(cell)){
      plot <- plot + ggraph::geom_node_text(aes(label=data, filter = depth == clonotype), color = label.color, size = label.size)
    } else {
      plot <- plot + ggraph::geom_node_text(aes(label=data, filter = depth == cell), color = label.color, size = label.size)
    }
  }

  if(label == "clonotype"){
    plot <- plot + ggraph::geom_node_text(aes(label=name, filter = depth == clonotype), color = label.color, size = label.size)
  }

  if(label == "specificity"){
    plot <- plot + ggraph::geom_node_text(aes(label=name, filter = depth == specificity), color = label.color, size = label.size)
  }

  return(plot)
}
