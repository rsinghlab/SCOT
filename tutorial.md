---
layout: template
---

## Tutorial
For many single-cell technologies, it is difficult to apply them simultaneously on the same single-cell. As a result, to get multiple views of a genome at the single-cell resolution, we need to integrate data from different -omic sequencing datasets. When these single-cell sequencing methods are applied to the same cell population or different populations that are expected to share some underlying biological manifold (e.g. common cell types), we expect there to be some cell-to-cell alignment to recover. SCOT performs alignment by using Gromov-Wasserstein optimal transport.  

<br> 
SCOT works in three steps:  
**1.** It first checks pairwise correlations between samples in each sequencing dataset and constructs k nearest-neighbor graphs for each to computing intra-domain distances while accounting for local geometry.  
**2.** It then optimizes the Gromov-Wasserstein optimal transport formulation to look for a probabilistic sample-wise correspondence matrix between the input datasets.  
**3.** Using this coupling matrix, it projects one dataset onto the other one using barycentric projection.  


![](assets/method_overview.png)

