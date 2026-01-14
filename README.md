# scMSGNN

**scMSGNN** (Multi-Scale Graph Neural Network for Single Cell Analysis) is an R package ported from the original Julia implementation of MSGNN. It leverages `torch` for R and the `Seurat` ecosystem to provide high-performance single-cell data analysis.

## Features

*   **Scalable Graph Learning**: Uses SIGN (Scalable Inception Graph Neural Networks) for efficient graph convolution without sampling.
*   **ZINB Denoising**: Built-in Zero-Inflated Negative Binomial autoencoder for robust denoising and imputation.
*   **Seurat Integration**: Seamlessly works with `Seurat` objects.
*   **Visualization**: Integration with `ggsc` and `ggtangle` for beautiful plots.

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("YuLab-SMU/scMSGNN")
```

## Usage

```r
library(scMSGNN)
library(Seurat)

# Load data
pbmc <- Read10X("path/to/data")
pbmc <- CreateSeuratObject(pbmc)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- RunPCA(pbmc)

# Run scMSGNN
pbmc <- RunscMSGNN(pbmc, 
                   k_neighbors = 20, 
                   sign_k = 3, 
                   epochs = 100, 
                   device = "cuda") # Use "cpu" if no GPU

# Visualize
plotCluster(pbmc, reduction = "msgnn")
```

## Requirements

*   R >= 4.0
*   Seurat >= 4.0
*   torch
