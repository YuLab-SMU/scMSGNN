#' Create Pathway Mask from GMT
#' 
#' Helper to create a binary mask for scMSGNN from a GMT file or list of pathways.
#' 
#' @param pathways List of pathways (e.g., read from GMT) where names are pathway names and values are gene vectors.
#' @param features Vector of gene names used in the model (must match order of input data).
#' @param sign_k Number of diffusion steps (to repeat the mask for diffused features).
#' @return A binary matrix [Pathways x (Genes * (sign_k + 1))]
#' @export
CreatePathwayMask <- function(pathways, features, sign_k) {
  n_genes <- length(features)
  n_pathways <- length(pathways)
  
  # Initialize mask for one scale: [Pathways x Genes]
  mask_base <- matrix(0, nrow = n_pathways, ncol = n_genes)
  rownames(mask_base) <- names(pathways)
  colnames(mask_base) <- features
  
  for (i in 1:n_pathways) {
    p_genes <- pathways[[i]]
    # Intersect with available features
    valid_genes <- intersect(p_genes, features)
    if (length(valid_genes) > 0) {
      mask_base[i, valid_genes] <- 1
    }
  }
  
  # Remove pathways with no genes? Optional.
  # For now keep them (all zeros) or user should filter.
  
  # Replicate mask for all scales (Original + sign_k diffused versions)
  # The input to model is [X, AX, A^2X ...]
  # We want the pathway node to connect to Gene_i in ALL scales.
  # So we repeat the mask horizontally.
  
  mask_full <- do.call(cbind, replicate(sign_k + 1, mask_base, simplify = FALSE))
  
  return(mask_full)
}
