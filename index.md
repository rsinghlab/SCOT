
![](assets/SCOT_logo.png)

**Note:** This GitHub page is a work in progress. It will contain information on the SCOT methodology, as well as demonstrated examples in various application cases.  

### Overview
SCOT is an unsupervised algorithm for cell-to-cell alignment of single-cell multi-omic data. <br>
For many single-cell technologies, it is difficult to apply them simultaneously on the same single-cell. As a result, to get multiple views of a genome at the single-cell resolution, we need to integrate data from different -omic sequencing datasets. When these single-cell sequencing methods are applied to the same cell population or different populations that are expected to share some underlying biological manifold (e.g. common cell types), we expect there to be some cell-to-cell alignment to recover. SCOT performs alignment by using Gromov-Wasserstein optimal transport.  

<br> 
SCOT works in three steps:  
**1.** It first checks pairwise correlations between samples in each sequencing dataset and constructs k nearest-neighbor graphs for each to computing intra-domain distances while accounting for local geometry.  
**2.** It then optimizes the Gromov-Wasserstein optimal transport formulation to look for a probabilistic sample-wise correspondence matrix between the input datasets.  
**3.** Using this coupling matrix, it projects one dataset onto the other one using barycentric projection.  

For more details, check out our [pre-print on bioRxiv](https://www.biorxiv.org/content/10.1101/2020.04.28.066787v2.full)

![](assets/method_overview.png)

### Example
