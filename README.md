---
output: github_document
---

# MuRD: Multi-Robust Cell-type Deconvolution of Bulk RNA Sequencing Data via Integrating Information from Single-cell Reference and External Expression Data
  
**MuRD** unifies multiple deconvolution schemes to infer cell type proportions from the target bulk RNA-seq data. Three unique features are embraced in this algorithm: first, **MuRD** is able to incorporate extra biological information from external data sources enables (e.g., scRNA_seq and other bulk RNA_seq data from independent studies); second, **MuRD** calibrates the reference-free algorithm by taking into account the proportion estimates from a reference-based approach; third, **MuRD** is robust to incorrect information from any of the provided data sources.

# Installation
You can install the released version of MuRD with:
```
#install devtools if necessary
install.packages("devtools")

#install the MuRD package
devtools::install_github("chencxxy28/MuRD")

#load
library(MuRD)
```

#Vignettes
Please visit [Tutorial](https://chencxxy28.github.io/MuRD/articles/MuRD.html).

