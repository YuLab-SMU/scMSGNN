#' @importFrom Matrix rowSums Diagonal t
#' @importFrom Seurat FindNeighbors
NULL

#' Normalize Adjacency Matrix
#'
#' Calculates the normalized Laplacian matrix: D^(-1/2) * (A + I) * D^(-1/2)
#'
#' @param adj A sparse adjacency matrix (dgCMatrix)
#' @return A normalized sparse matrix
#' @export
NormalizeDegree <- function(adj) {
  # 1. A + I (Self loops)
  # Check if diagonal is already 1 to avoid double adding
  diag(adj) <- 1
  
  # 2. Degree Matrix D
  D <- Matrix::rowSums(adj)
  
  # 3. D^(-1/2)
  D_inv_sqrt <- 1 / sqrt(D)
  D_inv_sqrt[is.infinite(D_inv_sqrt)] <- 0
  
  # 4. D^(-1/2) * A * D^(-1/2)
  # Efficient diagonal multiplication for sparse matrices
  D_mat <- Matrix::Diagonal(x = D_inv_sqrt)
  
  norm_adj <- D_mat %*% adj %*% D_mat
  return(norm_adj)
}

#' Build Graph from Seurat Object
#'
#' @param object Seurat object
#' @param k Number of neighbors
#' @param reduction Reduction to use (default: "pca")
#' @param dims Dimensions to use
#' @return A normalized adjacency matrix
#' @export
BuildGraph <- function(object, k = 20, reduction = "pca", dims = 1:30) {
  # Use Seurat's FindNeighbors to build KNN graph
  # Check if neighbors already exist for these parameters to save time? 
  # For now, just rebuild to be safe and consistent.
  
  object <- Seurat::FindNeighbors(object, 
                                  reduction = reduction, 
                                  dims = dims, 
                                  k.param = k, 
                                  verbose = FALSE)
  
  # Get the graph (usually stored as "RNA_snn" or similar)
  # Find the graph name
  graph_name <- paste0(Seurat::DefaultAssay(object), "_snn")
  if (!graph_name %in% names(object@graphs)) {
    # Fallback if name is different
    graph_name <- names(object@graphs)[1]
  }
  
  adj <- object@graphs[[graph_name]]
  
  # Normalize
  norm_adj <- NormalizeDegree(adj)
  
  return(norm_adj)
}
