#' SIGN Module
#' 
#' Multi-scale Graph Neural Network module
#' 
#' @export
nn_sign_module <- torch::nn_module(
  "SIGN",
  initialize = function(input_dim, hidden_dims, dropout = 0.0) {
    # input_dim: Dimension of concatenated features (Original + Diffused)
    
    layers <- list()
    
    in_d <- input_dim
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
