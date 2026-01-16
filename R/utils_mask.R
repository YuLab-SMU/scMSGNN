#' Create Pathway Mask from GMT
#' 
#' Helper to create a binary mask for scMSGNN from a GMT file or list of pathways.
#' 
#' @param gset.idx.list List of gene sets (e.g., read from GMT) where names are pathway/term names and values are gene vectors.
#' @param features Vector of gene names used in the model (must match order of input data).
#' @param sign_k Number of diffusion steps (to repeat the mask for diffused features).
#' @return A binary matrix [Pathways x (Genes * (sign_k + 1))]
#' @importFrom fastmatch %fin%
#' @export
CreatePathwayMask <- function(gset.idx.list, features, sign_k) {
  ## Initialize mask for one scale: [Pathways x Genes]

  ind <- lapply(gset.idx.list, function(i) which(features %fin% i))
  x <- Matrix::sparseMatrix(
         i = lapply(seq(length(ind)), function(i)rep(i, length(ind[[i]]))) |> unlist(),
         j = ind |> unlist(),
         x = 1,
         dims = c(length(gset.idx.list), length(features))
  )
  colnames(x) <- features
  rownames(x) <- names(gset.idx.list)

  # Remove pathways with no genes? Optional.
  # For now keep them (all zeros) or user should filter.
  
  # Replicate mask for all scales (Original + sign_k diffused versions)
  # The input to model is [X, AX, A^2X ...]
  # We want the pathway node to connect to Gene_i in ALL scales.
  # So we repeat the mask horizontally.
  
  x <- do.call(cbind, replicate(sign_k + 1, x, simplify = FALSE))
  
  return(x)
}
