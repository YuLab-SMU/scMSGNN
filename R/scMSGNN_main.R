#' Run scMSGNN
#'
#' Main entry point for running the MSGNN analysis.
#'
#' @param object Seurat object
#' @param features Features to use (default: VariableFeatures)
#' @param k_neighbors Number of neighbors for graph construction
#' @param sign_k Number of diffusion steps (SIGN k)
#' @param hidden_dims Hidden layer dimensions for MLP
#' @param epochs Number of training epochs
#' @param lr Learning rate
#' @param batch_size Batch size for training
#' @param device "cpu" or "cuda"
#' @param adj_matrix Optional: Custom adjacency matrix (sparse). If provided, skips graph building.
#' @param pathway_mask Optional: Binary mask matrix [Pathways x (Genes * (sign_k + 1))] for the first layer.
#' @importFrom torch optim_adam torch_index_select
#' @return Seurat object with "msgnn" reduction and denoised data
#' @export
RunscMSGNN <- function(object, 
                       features = NULL,
                       k_neighbors = 20,
                       sign_k = 2,
                       hidden_dims = c(128, 64), # Reduced for demo
                       epochs = 50,
                       lr = 0.001,
                       batch_size = 256,
                       device = "cpu",
                       reduction = 'pca',
                       reduction_dims = 30,
                       adj_matrix = NULL,
                       pathway_mask = NULL
                      ) {
  
  if (is.null(features)) {
    features <- Seurat::VariableFeatures(object)
  }
  if (length(features) == 0) {
    stop("No VariableFeatures found. Please run FindVariableFeatures first.")
  }

  if (torch::cuda_is_available() && device == "cuda") {
    device_obj <- torch::torch_device("cuda")
  } else {
    device_obj <- torch::torch_device("cpu")
  }
  
  # 1. Build Graph & Normalize
  if (is.null(adj_matrix)) {
    message("Building Graph...")
    # Assuming PCA is already run
    if (!"pca" %in% names(object@reductions)) {
      object <- Seurat::RunPCA(object, verbose = FALSE)
    }
    adj_matrix <- BuildGraph(object, k = k_neighbors, reduction = reduction, dims = seq(reduction_dims))
  } else {
    message("Using custom adjacency matrix...")
    # Ensure it's normalized or normalize it?
    # Let's assume user passes a raw adj and we normalize, or we check.
    # For safety, let's apply NormalizeDegree if it doesn't look normalized (hard to check cheaply).
    # Better: Assume user knows what they are doing OR apply normalization.
    # Let's apply NormalizeDegree which adds Self-Loop and normalizes.
    #norm_adj <- NormalizeDegree(adj_matrix)
  }

  adj_matrix <- ensure_tensor(adj_matrix, device = device_obj)
  adj_matrix <- normalize_adjacency(adj_matrix, device = device_obj)
  
  # 2. Prepare Data (Features)
  # X: Cells x Genes
  X <- Matrix::t(Seurat::GetAssayData(object, layer = "counts")[features, ])
  
  # Convert to Dense for Torch (Warning: Memory usage)
  # For very large datasets, we might need a custom dataset class
  X <- ensure_tensor(X, device = device_obj)
  
  message("Pre-computing Diffusion Features...")
  # 3. Pre-compute Diffusion: X, AX, A^2X ...
  # X is Cells x Genes. A is Cells x Cells.
  # Diffusion: A * X
  
  # Convert to Matrix for multiplication if not already
  # norm_adj is sparse
  
  feature_list <- list(X)
  current_X <- X
  
  for (i in 1:sign_k) {
    message(sprintf("  Diffusion step %d/%d", i, sign_k))
    current_X <- torch::torch_matmul(adj_matrix, current_X) 
    feature_list[[i + 1]] <- current_X
  }
  
  # Concatenate features: Cells x (Genes * (k+1))
  X_concat <- do.call(torch::torch_cat, list(feature_list, dim = 2))
  
  
  # 4. Initialize Model
  input_dim <- ncol(X_concat)
  output_dim <- ncol(X) # Predict original features
  
  # Handle Pathway Mask if provided
  mask_tensor <- NULL
  if (!is.null(pathway_mask)) {
    # Ensure mask is a matrix
    pathway_mask <- ensure_tensor(pathway_mask, device = device_obj)
    # Check dimensions
    if (ncol(pathway_mask) != input_dim) {
       stop(sprintf("Pathway mask columns (%d) must match input dimension (%d). Note: input_dim = n_genes * (sign_k + 1)", ncol(pathway_mask), input_dim))
    }
  }

  model <- nn_sign_module(input_dim, hidden_dims, pathway_mask = pathway_mask$to_dense())
  model$set_decoders(utils::tail(hidden_dims, 1), output_dim)
  model$to(device = device_obj)
  
  optimizer <- torch::optim_adam(model$parameters, lr = lr)
  loss_fn <- nn_zinb_loss()
  
  # 5. Training Loop
  message("Training SIGN Model...")
  
  num_samples <- nrow(X)
  num_batches <- ceiling(num_samples / batch_size)
  
  pb <- utils::txtProgressBar(min = 0, max = epochs, style = 3)
  
  for (epoch in 1:epochs) {
    total_loss <- 0
    
    # Shuffle indices
    indices <- sample(1:num_samples)
    
    for (b in 1:num_batches) {
      start_idx <- (b - 1) * batch_size + 1
      end_idx <- min(b * batch_size, num_samples)
      batch_idx <- indices[seq(start_idx, end_idx)]
      batch_idx <- torch::torch_tensor(batch_idx, dtype = torch::torch_long())
      
      batch_x <- torch_index_select(X_concat, 1, batch_idx)
      batch_target <- torch_index_select(X, 1, batch_idx) 
      
      optimizer$zero_grad()
      
      out <- model(batch_x)

      loss <- loss_fn(batch_target$to_dense(), out$mean, out$disp, out$pi)
      
      loss$backward()
      optimizer$step()

      total_loss <- total_loss + loss$item()
    }
    
    utils::setTxtProgressBar(pb, epoch)
  }
  close(pb)
  
  # 6. Inference (Get Embeddings and Denoised Data)
  model$eval()
  torch::with_no_grad({
    full_out <- model(X_concat$to(device = device_obj))
    embedding <- as.matrix(full_out$embedding$to(device = device_obj))
    denoised <- as.matrix(full_out$mean$to(device = device_obj))
  })
  
  rownames(embedding) <- colnames(object)
  colnames(embedding) <- paste0("MSGNN_", 1:ncol(embedding))
  
  rownames(denoised) <- colnames(object)
  colnames(denoised) <- features
  
  # 7. Store results in Seurat Object
  
  # Add Reduction
  object[["msgnn"]] <- Seurat::CreateDimReducObject(
    embeddings = embedding,
    key = "MSGNN_",
    assay = Seurat::DefaultAssay(object)
  )
  
  # Add Denoised Data (Optional: as a new Assay)
  # Transpose back to Genes x Cells
  object[["MSGNN"]] <- Seurat::CreateAssayObject(data = Matrix::t(denoised))
  
  message("Done.")
  return(object)
}
