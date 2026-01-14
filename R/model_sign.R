#' Masked Linear Layer
#' 
#' Applies a mask to the weights of a linear layer.
#' 
#' @export
nn_masked_linear <- torch::nn_module(
  "MaskedLinear",
  initialize = function(in_features, out_features, mask) {
    self$linear <- torch::nn_linear(in_features, out_features)
    # Register mask as a buffer (not a parameter to be trained)
    self$register_buffer("mask", mask)
  },
  
  forward = function(x) {
    # Apply mask to weights before linear transformation
    masked_weight <- self$linear$weight * self$mask
    torch::nn_functional_linear(x, masked_weight, self$linear$bias)
  }
)

#' SIGN Module
#' 
#' Multi-scale Graph Neural Network module
#' 
#' @export
nn_sign_module <- torch::nn_module(
  "SIGN",
  initialize = function(input_dim, hidden_dims, dropout = 0.0, pathway_mask = NULL) {
    # input_dim: Dimension of concatenated features (Original + Diffused)
    
    layers <- list()
    
    in_d <- input_dim
    
    # Check if we are using pathway mask for the first layer
    start_idx <- 1
    if (!is.null(pathway_mask)) {
      # First layer is Masked Linear
      # Mask shape: [out_features, in_features] (Torch Linear weight shape)
      # pathway_mask should be [n_pathways, n_genes * (k+1)]
      
      # Validate dimensions
      if (ncol(pathway_mask) != input_dim) {
        stop(sprintf("Mask input dim (%d) does not match data dim (%d)", ncol(pathway_mask), input_dim))
      }
      
      h_d <- nrow(pathway_mask) # Number of pathways
      
      layers[[length(layers) + 1]] <- nn_masked_linear(in_d, h_d, pathway_mask)
      layers[[length(layers) + 1]] <- torch::nn_batch_norm1d(h_d)
      layers[[length(layers) + 1]] <- torch::nn_relu()
      if (dropout > 0) {
        layers[[length(layers) + 1]] <- torch::nn_dropout(dropout)
      }
      
      in_d <- h_d
      # If hidden_dims provided, continue building subsequent layers
      # Note: The first hidden dim in `hidden_dims` is ignored if mask is used, 
      # or we expect user to align them. 
      # Let's assume `hidden_dims` defines layers AFTER the pathway layer if mask is present.
    }
    
    for (h_d in hidden_dims) {
      layers[[length(layers) + 1]] <- torch::nn_linear(in_d, h_d)
      layers[[length(layers) + 1]] <- torch::nn_batch_norm1d(h_d)
      layers[[length(layers) + 1]] <- torch::nn_relu()
      if (dropout > 0) {
        layers[[length(layers) + 1]] <- torch::nn_dropout(dropout)
      }
      in_d <- h_d
    }
    
    self$encoder <- torch::nn_sequential(layers)
    
    # ZINB Decoders
    # We need to predict Mean, Disp, Pi for original dimension
    # But usually AutoEncoders reconstruct the input.
    # In SIGN, the input to MLP is concatenated diffused features.
    # The output target is usually the original expression (or denoised version).
    # Let's assume the output dimension is 'input_dim / (k+1)' (original features) 
    # OR we pass output_dim explicitly.
    # Here we assume the last hidden layer maps to 3 heads: Mean, Disp, Pi
  },
  
  set_decoders = function(last_hidden_dim, output_dim) {
    self$dec_mean <- torch::nn_sequential(
      torch::nn_linear(last_hidden_dim, output_dim),
      torch::nn_softplus() # Mean must be positive
    )
    self$dec_disp <- torch::nn_sequential(
      torch::nn_linear(last_hidden_dim, output_dim),
      torch::nn_softplus() # Dispersion must be positive
    )
    self$dec_pi <- torch::nn_sequential(
      torch::nn_linear(last_hidden_dim, output_dim),
      torch::nn_sigmoid() # Probability [0, 1]
    )
  },
  
  forward = function(x) {
    # x: [Batch, Input_Dim] (Concatenated features)
    enc <- self$encoder(x)
    
    mean <- self$dec_mean(enc)
    disp <- self$dec_disp(enc)
    pi <- self$dec_pi(enc)
    
    return(list(mean = mean, disp = disp, pi = pi, embedding = enc))
  }
)
