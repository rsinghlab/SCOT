---
layout: template
---

![](assets/SCOT_logo.png)

### Welcome to the documentation page for SCOT!
**Note:** This website is currently under construction; it will be complete by the end of Dec 2020.<br>

SCOT is an unsupervised algorithm for performing cell-to-cell alignment of single-cell multi-omic datasets. <br>
We provide tutorials and examples on this website. For details on methodology and experimental results, please check out [our pre-print on bioRxiv](https://www.biorxiv.org/content/10.1101/2020.04.28.066787v2.full)

### Overview
Availability of various sequencing technologies allows us to capture different properties of the genome at the single-cell resolution, such as transcriptomic and  epigenomic states. This is useful to get a holistic view of the single-cell genome and help study the molecular mechanisms behind gene regulation or cellular heterogeneity from different genomic aspects. <br>

However, with very few exceptions of co-assaying technologies, applying different sequencing assays on the same single-cell is not possible. For example, for some combinations of the sequencing assays that need to access to the same part of the genome (e.g. sci-Hi-C and sci-ATAC-seq), there are currently no co-assaying technology available. In addition, one cannot sequentially apply these methods on the same cell either because a cell is destroyed after any sequencing procedure.<br>

As a result, to get a holistic view of the genome, we tend to divide a cell population of interest into different aliquots, and then apply a different sequencing assay on each one. Then, to perform a joint multi-omic analysis, we need to integrate these datasets. As long as the measurements have been taken from the cell populations that are expected to have some shared underlying biological manifold (e.g. cell populations that contain the same cell types and states or aliquots that come from the same population, where cells share geneology), we expect to have sample-wise correspondences to recover across measurements and be able to align these datasets for an integrated view.<br>

![](assets/problem.png)

However, this is a challenging problem that requires an unsupervised computational approach because the process yields disparate datasets. Since we measure different properties of the genome on different cells, we have no 1-to-1 correspondence information between either the cells or the features. <br>

**SCOT** is an unsupervised alignment tool that yields correspondence probabilities between cells from different -omics datasets. Unlike many other integration methods, it performs alignment without requiring any correspondence information *a priori* and aligns datasets based on the correspondence probabilities it recovers.

#### SCOT works in three steps:  
**1.** It first checks pairwise correlations between samples in each sequencing dataset and constructs k nearest-neighbor (kNN) graphs for each. Based on the shortest distances on these kNN graphs, it computes intra-domain distances.
**2.** It then optimizes the Gromov-Wasserstein optimal transport formulation to look for a probabilistic correspondence matrix between the samples of the input datasets that will minimize the Euclidean distances between the intra-domain distance matrices. 
**3.** Using this correspondence matrix, it projects one dataset onto the other via barycentric projection.  

![](assets/method_overview.png)

### Why should I prefer SCOT to use on my dataset?
While there are other alignment tools available, SCOT brings a few advantages that are important in real-world settings:<br>
**1.** Most alignment methods are developed for batch integration of single-cell RNA-seq datasets (e.g. scAlign, MNN, Seurat, Harmony) and are [shown to perform poorly on multi-omic alignment tasks](https://academic.oup.com/bioinformatics/article/36/Supplement_1/i48/5870490#206061395), which is fundamentally a different problem. SCOT is specifically designed for and [tested on multi-omic integration tasks](https://www.biorxiv.org/content/10.1101/2020.04.28.066787v2).<br>
**2.** Many alignment tools (e.g. MAGAN, Seurat) require anchor points or feature-wise correspondence information to perform alignment, which one does not likely to have upon sequencing separate cell populations. SCOT is an **unsupervised** tool and does not require any correspondence information to be known *a priori* in order to perform alignment.<br>
**3.** Other currently available unsupervised multi-omic alignment tools (UnionCom, MMD-MA, Pamona) require users to perform hyperparameter optimization in order to yield high quality alignments. Without any validation data on correspondences, it is difficult to perform hyperparameter tuning. SCOT provides [a procedure to guide hyperparameter selection in fully unsupervised settings](https://rsinghlab.github.io/SCOT/unsupervised/). <br>
**4.** In comparison to the other unsupervised mutli-omic alignment tools, SCOT is [computationally scalable to large datasets](https://www.biorxiv.org/content/10.1101/2020.04.28.066787v2.full.pdf).<br>


If you would like to use SCOT for your datasets, take a look at [the tutorial page](rsinghlab.github.io/SCOT/tutorial), as well as [the examples]() we provide on this documentation site. If you are having trouble getting started or have questions about SCOT, please do not hesitate to [**contact us**](rsinghlab.github.io/SCOT/contact)!





