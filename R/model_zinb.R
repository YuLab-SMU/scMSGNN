#' ZINB Loss Module
#' 
#' Implements the Zero-Inflated Negative Binomial Loss
#' 
#' @export
nn_zinb_loss <- torch::nn_module(
  "ZINBLoss",
  initialize = function() {
    # No parameters to initialize for the loss itself
  },
  
  forward = function(x, mean, disp, pi, scale_factor = 1.0) {
    # x: Target data (Observed counts)
    # mean: Predicted mean (mu)
    # disp: Predicted dispersion (theta)
    # pi: Predicted dropout probability (pi)
    
    eps <- 1e-7
    
    # Scale mean by library size if provided (or handled outside)
    # Here we assume 'mean' is already the predicted count parameter
    
    # Case 1: x = 0
    # Prob = pi + (1-pi) * NB(0)
    # NB(0) = (theta / (theta + mu))^theta
    
    t1 <- torch::torch_pow(disp / (disp + mean + eps), disp)
    zero_case <- -torch::torch_log(pi + ((1.0 - pi) * t1) + eps)
    
    # Case 2: x > 0
    # Prob = (1-pi) * NB(x)
    # LogProb = log(1-pi) + log(NB(x))
    # log(NB(x)) = lgamma(x+theta) - lgamma(theta) - lgamma(x+1) + theta*log(theta/(theta+mu)) + x*log(mu/(theta+mu))
    
    t2 <- torch::torch_lgamma(disp + eps) + torch::torch_lgamma(x + 1.0) - torch::torch_lgamma(x + disp + eps)
    
    # Correction: The formula in Julia code:
    # lg(Θ+eps) + lg(X+1f0) - lg(X+Θ+eps)  <-- This seems to be -log(Combination)? 
    # Standard NB LogProb:
    # lgamma(x+theta) - lgamma(x+1) - lgamma(theta) + theta * log(theta) + x * log(mu) - (x+theta) * log(theta+mu)
    
    # Let's follow the standard PyTorch ZINB implementation or the Julia one carefully.
    # Julia:
    # t1 = lg(Θ) + lg(X+1) - lg(X+Θ) 
    # t2 = (Θ+X) * log(1 + (M/Θ))  <-- log((Θ+M)/Θ) = - log(Θ/(Θ+M))
    # t3 = X * (log(Θ) - log(M))   <-- X * log(Θ/M)
    # This looks like a specific rearrangement.
    
    # Let's use the standard formulation which is numerically stable:
    # log_prob = lgamma(x + theta) - lgamma(theta) - lgamma(x + 1)
    #          + theta * (log(theta) - log(theta + mu))
    #          + x * (log(mu) - log(theta + mu))
    
    nb_case <- torch::torch_lgamma(x + disp + eps) - torch::torch_lgamma(disp + eps) - torch::torch_lgamma(x + 1.0) +
               disp * (torch::torch_log(disp + eps) - torch::torch_log(disp + mean + eps)) +
               x * (torch::torch_log(mean + eps) - torch::torch_log(disp + mean + eps))
    
    non_zero_case <- - (torch::torch_log(1.0 - pi + eps) + nb_case)
    
    # Combine
    zero_mask <- (x < eps)
    loss <- torch::torch_where(zero_mask, zero_case, non_zero_case)
    
    return(torch::torch_mean(loss) * scale_factor)
  }
)
