---
output: github_document
---

# MultiRD: Multi-Robust Cell-type Deconvolution of Bulk RNA Sequencing Data via Integrating Information from Single-cell Reference and External Expression Data
  
**MultiRD** unifies multiple deconvolution schemes to infer cell type proportions from the target bulk RNA-seq data. Three unique features are embraced in this algorithm: first, **MultiRD** is able to incorporate extra biological information from external data sources enables (e.g., scRNA_seq and other bulk RNA_seq data from independent studies); second, **MultiRD** calibrates the reference-free algorithm by taking into account the proportion estimates from a reference-based approach; third, **MultiRD** is robust to incorrect information from any of the provided data sources.

# Installation
You can install the released version of MultiRD with:
```
#install devtools if necessary
install.packages("devtools")

#install the MultiRD package
devtools::install_github("chencxxy28/MultiRD")

#load
library(MultiRD)
```

If [_Biobase_](https://bioconductor.org/packages/release/bioc/html/Biobase.html) package is not available, please install it first before installation of **MultiRD**
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biobase")
```

# Vignettes
Please visit [Tutorial](https://chencxxy28.github.io/MultiRD/articles/MultiRD.html).

