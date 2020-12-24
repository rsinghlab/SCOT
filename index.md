---
layout: template
---

![](assets/SCOT_logo.png)

### Welcome to the documentation page for SCOT! <br>

SCOT is an unsupervised algorithm for performing cell-to-cell alignment on single-cell multi-omic datasets. <br>
We provide tutorials and examples on this website. For details on methodology and experimental results, please check out [our pre-print on bioRxiv](https://www.biorxiv.org/content/10.1101/2020.04.28.066787v2.full)

### Overview
Thanks to the availability of various sequencing technologies, we are able to capture different properties of the genome at the single-cell resolution, such as transcriptomic and different epigenomic states. This is useful to get a holistic view of the single-cell genome and perhaps study the molecular mechanisms behind gene regulation, or cellular heterogeneity from different genomic aspects. <br>

However, with very few exceptions, applying different sequencing assays on the same single-cell is difficult. For some combinations of the sequencing assays that need to access to the same part of the genome (e.g. sci-Hi-C and sci-ATAC-seq), this is currently not possible at all. In addition, one cannot sequentially apply these methods on the same cell either because a cell is destroyed after any sequencing procedure.<br>


As a result, to get a holistic view of the genome, we tend to divide a cell population of interest into different aliquots, and then apply a different sequencing assay on each one. Then, to perform a joint multi-omic analysis, we need to integrate these datasets. 

![](assets/problem2.png)

When these single-cell sequencing methods are applied to the same cell population or different populations that are expected to share some underlying biological manifold (e.g. common cell types), we expect there to be some cell-to-cell alignment to recover.  

However, this is a challenging problem because the process yields disparate datasets.



SCOT is an unsupervised alignment algorithm that recovers probabilistic correspondences between the cells of different -omic datasets. It is appropriate to apply SCOT for multi-omic alignment as long as you have a reason to believe there ... 

SCOT yields probabilistic cell-to-cell correspondences and performs alignment by using Gromov-Wasserstein optimal transport.  


The image above shows an example of 

that requires unsupervised alignment algorithms 

First of all, we have no correspondence information between cells , and often, we don’t have that for features either since we just measured different properties of genome in different cells.  So, this process yields disparate datasets with no information on alignment a priori. As a result, we need unsupervised algorithms that will align them without relying on any correspondence information.
 <br>



#### SCOT works in three steps:  
**1.** It first checks pairwise correlations between samples in each sequencing dataset and constructs k nearest-neighbor (kNN) graphs for each. Based on the shortest distances on these kNN graphs, it computes intra-domain distances.
**2.** It then optimizes the Gromov-Wasserstein optimal transport formulation to look for a probabilistic sample-wise correspondence matrix between the input datasets.  
**3.** Using this correspondence matrix, it projects one dataset onto the other one via barycentric projection.  

![](assets/method_overview.png)

### Why should I use SCOT on my dataset?
While there are other alignment tools available, SCOT brings a few advantages that we feel are important in real-world settings:
**1.** Most alignment methods are developed for batch integration of single-cell RNA-seq datasets (scAlign, MNN, Seurat) and are [shown to perform poorly on multi-omic alignment tasks](). SCOT is 
**2.** Some of the alignment tools require anchor points or feature-wise correspondence information to perform alignment, which ... ... SCOT does not require any correspondence information to be known *a priori* in order to perform alignment.
**3.** Most of the unsupervised multi-omic alignment tools require hyperparameter optimization in order to ... Although 
**4.** 

SCOT yields alignments 

If you would like to use SCOT for your datasets but are having trouble getting started after following the [tutorial](rsinghlab.github.io/SCOT/tutorial) or the [examples](rsinghlab.github.io/SCOT/examples), or have questions about it, please do not hesitate to [**contact us**](rsinghlab.github.io/SCOT/contact)!





