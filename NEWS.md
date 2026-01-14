# scMSGNN 0.2.0

## New Features (Interpretable AI)

*   **Pathway Masking Support**: Added `pathway_mask` argument to `RunscMSGNN()`. This allows users to constrain the neural network architecture based on biological pathways (e.g., KEGG, Reactome), making the deep learning model interpretable.
*   **Custom Graph Topology**: Added `adj_matrix` argument to `RunscMSGNN()`, enabling the integration of external biological knowledge (e.g., protein-protein interaction networks or ligand-receptor interactions) into the graph convolution process.
*   **New Helper Functions**:
    *   `CreatePathwayMask()`: A utility to generate binary mask matrices from gene sets for the neural network.
    *   `nn_masked_linear()`: A new PyTorch module that implements a linear layer with frozen non-pathway connections.

## Improvements

*   Enhanced `RunscMSGNN` to automatically handle device selection and mask tensor conversion.

# scMSGNN 0.1.0

## Initial Release

*   **Port from Julia**: This package is a complete port of the original `MSGNN` algorithm (implemented in Julia) to R, leveraging the `torch` package for backend computations and `Seurat` for data handling.
*   **Scalable Graph Learning**: Implemented the SIGN (Scalable Inception Graph Neural Networks) architecture, allowing for efficient training on large-scale single-cell datasets without the need for graph sampling.
*   **ZINB Denoising**: Integrated a Zero-Inflated Negative Binomial (ZINB) autoencoder to robustly model sparse single-cell RNA-seq data and recover true expression signals.
*   **Seurat Integration**: Seamless workflow integration with `Seurat` objects, adding `msgnn` reduction and denoised assays directly to the object.
