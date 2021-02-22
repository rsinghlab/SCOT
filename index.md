---
layout: template
---

![](assets/SCOT_logo.png)

### Welcome to the documentation page for SCOT!
SCOT is an unsupervised algorithm for performing cell-to-cell alignment of single-cell multi-omic datasets. <br>
We provide tutorials and examples on this website. For details on methodology and experimental results, please check out [our pre-print on bioRxiv](https://www.biorxiv.org/content/10.1101/2020.04.28.066787v2.full)

### Overview
Availability of various sequencing technologies allows us to capture different properties of the genome at the single-cell resolution, such as transcriptomic and  epigenomic states, and get a holistic view of the single-cell genome. However, with very few exceptions of co-assaying technologies, applying different sequencing assays on the same single-cell is not possible. One way scientist can get around this is by dividing a cell population of interest into different aliquots, and then applying a different sequencing assay on each one. Then, to perform a joint multi-omic analysis, they need to integrate these datasets. As long as the measurements have been taken from the cell populations that are expected to have some shared underlying biological manifold (e.g. cell populations that contain the same cell types/states or aliquots that come from the same population, where cells share geneology), we expect to have probabilistic sample-wise correspondences to recover across measurements and be able to align these datasets for an integrated view.<br>

This problem requires unsupervised computational methods because the process yields disparate datasets: Since we measure different properties of the genome on different cells, we have no 1-to-1 correspondence information between either the cells or the features. <br>
![](assets/problem.png)

**SCOT** is an **unsupervised alignment tool** that yields correspondence probabilities between cells from different -omics datasets. Unlike many other integration methods, it runs **without requiring any correspondence information** *a priori* and aligns datasets based on the correspondence probabilities it recovers.

#### SCOT works in three steps:  
**1.** It first checks pairwise correlations between samples in each sequencing dataset and constructs k nearest-neighbor (kNN) graphs for each. Based on the shortest distances on these kNN graphs, it computes intra-domain distances.
**2.** It then optimizes the Gromov-Wasserstein optimal transport formulation to look for a probabilistic correspondence matrix between the samples of the input datasets that will minimize the Euclidean distances between the intra-domain distance matrices. 
**3.** Using this correspondence matrix, it projects one dataset onto the other via barycentric projection.  

![](assets/method_overview.png)

### Why should I prefer SCOT to use on my dataset?
While there are other alignment tools available, SCOT brings a few advantages that are important in real-world settings:<br>
**1. High quality multi-omic alignment** Most alignment methods are developed for batch integration of single-cell RNA-seq datasets (e.g. scAlign, MNN, Seurat, Harmony) and are [shown to perform poorly on multi-omic alignment tasks](https://academic.oup.com/bioinformatics/article/36/Supplement_1/i48/5870490#206061395), which is fundamentally a different problem. **SCOT is specifically designed for and [tested on multi-omic integration tasks](https://www.biorxiv.org/content/10.1101/2020.04.28.066787v2).**<br>
**2. Unsupervised alignment** Many alignment tools (e.g. MAGAN, Seurat) require anchor points or feature-wise correspondence information to perform alignment, which one does not likely to have upon sequencing separate cell populations. **SCOT is an unsupervised tool and does not require any correspondence information to be known *a priori* in order to perform alignment.**<br>
**3. Approximately self-tuning hyperparameters** Other currently available unsupervised multi-omic alignment tools (UnionCom, MMD-MA, Pamona) require users to perform hyperparameter optimization in order to yield high quality alignments. Without any validation data on correspondences, it is difficult to perform hyperparameter tuning. **SCOT provides [a procedure to approximately self-tune hyperparameters in fully unsupervised settings](https://rsinghlab.github.io/SCOT/unsupervised/).** <br>
**4. Handling cell type imbalance** Through an extension with unbalanced optimal transport, **SCOT is able to [handle cell type imbalance between multi-omic assays](https://rsinghlab.github.io/SCOT/examples_unbalanced).**
**5. Computational scalability** In comparison to the other unsupervised mutli-omic alignment tools, SCOT is [computationally scalable to large datasets](https://www.biorxiv.org/content/10.1101/2020.04.28.066787v2.full.pdf).<br>

#### Contact and Citation
If you would like to use SCOT for your datasets, take a look at [the tutorial page](rsinghlab.github.io/SCOT/tutorial), as well as [the examples]() we provide on this documentation site. If you are having trouble getting started or have questions about SCOT, please do not hesitate to [**contact us**](rsinghlab.github.io/SCOT/contact)!<br>
If you use SCOT for your work, you can cite our pre-print as below:<br>
**BibTex Citation:**  
```
@article {SCOT2020,  
	author = {Demetci, Pinar and Santorella, Rebecca and Sandstede, Bj{\"o}rn and Noble, William Stafford and Singh, Ritambhara},  
	title = {Gromov-Wasserstein optimal transport to align single-cell multi-omics data},  
	elocation-id = {2020.04.28.066787},  
	year = {2020},  
	doi = {10.1101/2020.04.28.066787},  
	publisher = {Cold Spring Harbor Laboratory},  
	URL = {https://www.biorxiv.org/content/early/2020/11/11/2020.04.28.066787},  
	eprint = {https://www.biorxiv.org/content/early/2020/11/11/2020.04.28.066787.full.pdf},  
	journal = {bioRxiv}. 
}
```



