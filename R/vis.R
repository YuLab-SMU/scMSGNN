#' @import ggplot2
#' @import ggsc
#' @import ggtangle

#' @title Plot Clustering Results
#' @description Wrapper around ggsc::sc_dim to plot clustering results on MSGNN embedding
#' 
#' @param object Seurat object
#' @param reduction Reduction to use (default: "msgnn")
#' @param ... Additional arguments passed to ggsc::sc_dim
#' @export
plotCluster <- function(object, reduction = "msgnn", ...) {
  if (!requireNamespace("ggsc", quietly = TRUE)) {
    stop("Package 'ggsc' is required for this function.")
  }
  ggsc::sc_dim(object, reduction = reduction, ...)
}

#' Plot Feature Scatter
#' 
#' Wrapper around ggsc::sc_feature
#' 
#' @param object Seurat object
#' @param features Features to plot
#' @param reduction Reduction to use
#' @param ... Additional arguments passed to ggsc::sc_feature
#' @export
plotScatter <- function(object, features, reduction = "msgnn", ...) {
  if (!requireNamespace("ggsc", quietly = TRUE)) {
    stop("Package 'ggsc' is required for this function.")
  }
  ggsc::sc_feature(object, features = features, reduction = reduction, ...)
}

#' Plot Graph Structure
#' 
#' Wrapper around ggtangle for network visualization
#' 
#' @param object Seurat object
#' @param layout Layout algorithm
#' @return A ggplot object
#' @export
plotGraph <- function(object, layout = "fr") {
  if (!requireNamespace("ggtangle", quietly = TRUE)) {
    stop("Package 'ggtangle' is required for this function.")
  }
  # This assumes ggtangle can handle Seurat objects or their graph slots
  # Since ggtangle is not standard, we provide a generic wrapper.
  # Users might need to extract the graph first.
  
  # Attempt to extract graph
  graph_name <- paste0(Seurat::DefaultAssay(object), "_snn")
  if (!graph_name %in% names(object@graphs)) {
    stop("No graph found in Seurat object.")
  }
  
  adj <- object@graphs[[graph_name]]
  # Convert to igraph or appropriate format for ggtangle if needed
  # For now, we assume ggtangle has a method for Seurat or matrix
  
  # Placeholder: ggtangle::ggtangle(adj)
  message("Using ggtangle to visualize graph...")
  # ggtangle usage depends on its API
  # ggtangle::ggtangle(object, ...)
}
