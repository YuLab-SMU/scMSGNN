#' Check if object is a torch tensor
#'
#' @param x Object to check
#' @return Logical
#' @export
is_torch_tensor <- function(x) {

  inherits(x, "torch_tensor")
}

#' Convert to torch tensor if needed
#'
#' @param x Matrix or tensor
#' @param device Device to use ("cpu" or "cuda")
#' @param ... additional parameters
#' @importFrom torch torch_float
#' @return torch_tensor
#' @export
ensure_tensor <- function(x, device = 'cpu', ...) {
  UseMethod("ensure_tensor")
}

#' @method ensure_tensor torch_tensor
#' @export
ensure_tensor.torch_tensor <- function(x, device = 'cpu', ...){
   return(x$to(device = device))
}

#' @method ensure_tensor matrix
#' @export
ensure_tensor.matrix <- function(x, device = 'cpu', ...){
   torch_tensor(x, dtype = torch_float())$to(device = device)
}

#' @method ensure_tensor data.frame
#' @export
ensure_tensor.data.frame <- function(x, device = 'cpu', ...){
   torch_tensor(as.matrix(x), dtype = torch_float())$to(device = device)
}

#' @importFrom torch torch_tensor torch_sparse_coo_tensor torch_long
#' @method ensure_tensor sparseMatrix
#' @export
ensure_tensor.sparseMatrix <- function(x, device = 'cpu', ...){
  x_coo <- as(x, "TsparseMatrix")
  if (Matrix::isSymmetric(x_coo)){
    indi <- c(x_coo@i, x_coo@j) + 1L
    indj <- c(x_coo@j, x_coo@i) + 1L
    indices <- torch_tensor(rbind(indi, unlist(indj)), dtype = torch_long())
    values <- rep(x_coo@x, 2)
  }else{
    indices <- torch_tensor(rbind(x_coo@i + 1L, x_coo@j + 1L), dtype = torch_long())
    values <- torch_tensor(x_coo@x, dtype = torch_float())
  }

  shape <- x_coo@Dim

  sparse_ts <- torch_sparse_coo_tensor(indices, values, size = shape)

  return(sparse_ts$to(device = device))

}



#' Normalize adjacency matrix for GCN
#'
#' @description Computes D^(-0.5) * A * D^(-0.5)
#' @param adj Adjacency matrix (tensor or matrix)
#' @param add_self_loops Whether to add self-loops
#' @param device Device to use
#' @importFrom torch torch_eye torch_sum torch_pow torch_where torch_matmul
#' @importFrom torch torch_diag torch_isinf torch_zeros_like torch_ones
#' @return Normalized adjacency tensor
#' @export
normalize_adjacency <- function(adj, add_self_loops = TRUE, device = "cpu") {
  adj <- ensure_tensor(adj, device)
  n <- adj$size(1)

  if (!adj$is_sparse()){
    if (add_self_loops) {
      eye <- torch_eye(n, device = device)
      adj <- adj + eye
    }

    # Degree matrix
    deg <- torch_sum(adj, dim = 2)
    deg_inv_sqrt <- torch_pow(deg + 1e-10, -0.5)
    deg_inv_sqrt <- torch_where(
      torch_isinf(deg_inv_sqrt),
      torch_zeros_like(deg_inv_sqrt),
      deg_inv_sqrt
    )
    d_mat <- torch_diag(deg_inv_sqrt)

    # Symmetric normalization
    return(torch_matmul(torch_matmul(d_mat, adj), d_mat))
  }

  return(.normalize_adjacency_sparse(adj, n, add_self_loops))
}



.normalize_adjacency_sparse <- function(adj, n, add_self_loops = TRUE){
  adj <- adj$coalesce()
  indices <- adj$indices()
  values <- adj$values()

  if (add_self_loops) {
    loop_indices <- torch_tensor(rbind(seq(n), seq(n)), dtype = torch_long(), device = adj$device)
    loop_values <- torch_ones(n, dtype = values$dtype, device = adj$device)
    loop_adj <- torch_sparse_coo_tensor(loop_indices, loop_values, c(n, n), device=adj$device)$coalesce()

    adj <- adj + loop_adj
    indices <- adj$indices()
    values <- adj$values()
  }
  # Calculate Degree: Sum rows
  # Sparse sum returns a dense tensor
  ones_vec <- torch_ones(n, 1, dtype = values$dtype, device = adj$device)
  deg <- torch_matmul(adj, ones_vec)$squeeze(2)
  deg_inv_sqrt <- torch_pow(deg + 1e-10, -0.5)
  deg_inv_sqrt <- torch_where(
    torch_isinf(deg_inv_sqrt),
    torch_zeros_like(deg_inv_sqrt),
    deg_inv_sqrt
  )
  row <- indices[1, ]
  col <- indices[2, ]

  d_row <- deg_inv_sqrt$index_select(1, row + 1L)
  d_col <- deg_inv_sqrt$index_select(1, col + 1L)

  norm_values <- values * d_row * d_col

  torch_sparse_coo_tensor(indices + 1L, norm_values, c(n, n))$coalesce()

}

