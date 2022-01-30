
![](assets/SCOT_logo.png)

# Single-Cell alignment using Optimal Transport (SCOT)
* **Note:** We are currently updating this repository with compertmantalized versions of the algorithm due to the new development of v.2.0.

SCOT is a Python tool for performing unsupervised alignment of single-cell multi-omics datasets. Its methodology is detailed in the following two papers:
- SCOT v.1.0: [Gromov-Wasserstein based optimal transport for aligning single-cell multi-omics data](https://www.biorxiv.org/content/10.1101/2020.04.28.066787v2)
- SCOT v.2.0: [Unsupervised integration of single-cell multi-omics datasets with disparities in cell-type representation](https://www.biorxiv.org/content/10.1101/2021.11.09.467903v1)

For full documentation, please visit https://rsinghlab.github.io/SCOT/ (currently being updated).  

## SCOT v.1.0
Unsupervised single-cell multi-omic integration with Gromov-Wasserstein optimal transport & a self-tuning heuristic for hyperparameter selection.

***THIS ARCHIVE CONTAINS SCOT v.1.0***
Usage: All dependencies are recorded in `requirements.txt`. You can install them together with `pip install requirements.txt`.
Jupyter notebooks to replicate the results from the manuscript are in the folder `/replication`. These also give examples for how to use SCOT. Scripts in `/examples` contain sample scripts for unsupervised and supervised hyperparameter selection. 
E-mail: `pinar_demetci@brown.edu`, `pinardemetci@gmail.com`, `rebecca_santorella@brown.edu` or `ritambhara@brown.edu` if you have any questions.

Basic use:
```{python} 
# Given two numpy matrices, domain1 and domain2, where the rows are cells and columns are different genomic features:
scot= SCOT(domain1, domain2)
aligned_domain1, aligned_domain2 = scot.align(k=50, e=1e-3)

#If you can't pick the parameters k and e, you can try out our unsupervised self-tuning heuristic by running:
scot= SCOT(domain1, domain2)
aligned_domain1, aligned_domain2 = scot.align(selfTune=True)
```
**Required parameters for the `align` method:**
- *k:* Number of neighbors to be used when constructing kNN graphs. Default= min(min(n_1, n_2), 50), where n_i, for i=1,2 corresponds to the number of samples in the i^th domain.
- *e:* Regularization constant for the entropic regularization term in entropic Gromov-Wasserstein optimal transport formulation. Default= 1e-3 
   
**Optional parameters:**
- *normalize=* Determines whether to normalize input data ahead of alignment. True or False (boolean parameter). Default = True.
- *norm=* Determines what sort of normalization to run, "l2", "l1", "max", "zscore". Default="l2" 
- *mode:* "connectivity" or "distance". Determines whether to use a connectivity graph (adjacency matrix of 1s/0s based on whether nodes are connected) or a distance graph (adjacency matrix entries weighted by distances between nodes). Default="connectivity"  
- *metric:* Sets the metric to use while constructing nearest neighbor graphs. some possible choices are "correlation", "minkowski".  "correlation" is Pearson's correlation and "minkowski" is equivalent to Euclidean distance in its default form (). Default= "correlation". 
- *verbose:* Prints loss while optimizing the optimal transport formulation. Default=True
- *XontoY:* Determines the direction of barycentric projection. True or False (boolean parameter). If True, projects domain1 onto domain2. If False, projects domain2 onto domain1. Default=True.

***Note:*** If you want to specify the marginal distributions of the input domains and not use uniform distribution, please set the attributes p and q to the distributions of your choice (for domain 1, and 2, respectively) after initializing a SCOT class instance and before running alignment and set init_marginals=False in .align() parameters

## SCOT v.1.1
A naive extension to multi-modal alignment, where the first dataset in the input as treated as the anchor to align on. 

## SCOT v.2.0
A few extensions:
1) Alignment with the unbalanced Gromov-Wasserstein optimal transport formulation to handle cell-type representation disparities (Sejourne et al, 2020)
2) Multi-modal alignment by picking the anchor domain based on imputation potential of domain-specific nearest neighbor graphs
3) Different choices for joint embedding/projection

### Citation:
We are excited to see any extentions and improvements our work! If you are using code from this repository, please kindly cite our work: 

**For SCOT v.1.0:** <br>
Demetci, P. Santorella, R. Sandstede, B., Noble, W. S., Singh, R. 2020. Gromov-Wasserstein based optimal transport for aligning single-cell multi-omics data. bioRxiv. 2020.04.28.066787; doi: https://doi.org/10.1101/2020.04.28.066787
**BibTex Citation:**  
```
@article {Demetci2020.SCOT,  
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

**For SCOT v.2.0:** <br>
Demetci, P. Santorella, R. Sandstede, B., Noble, W. S., Singh, R. 2021. Unsupervised integration of single-cell multi-omics datasets with disparities in cell-type representation. bioRxiv. 2021.11.09.467903; doi: https://doi.org/10.1101/2021.11.09.467903
**BibTex Citation:**  
```

@article{Demetci2021.SCOTv2,
	author = {Demetci, Pinar and Santorella, Rebecca and Sandstede, Bj{\"o}rn and Singh, Ritambhara},
	doi = {10.1101/2021.11.09.467903},
	elocation-id = {2021.11.09.467903},
	eprint = {https://www.biorxiv.org/content/early/2021/11/11/2021.11.09.467903.full.pdf},
	journal = {bioRxiv},
	publisher = {Cold Spring Harbor Laboratory},
	title = {Unsupervised integration of single-cell multi-omics datasets with disparities in cell-type representation},
	url = {https://www.biorxiv.org/content/early/2021/11/11/2021.11.09.467903},
	year = {2021},
	Bdsk-Url-1 = {https://www.biorxiv.org/content/early/2021/11/11/2021.11.09.467903},
	Bdsk-Url-2 = {https://doi.org/10.1101/2021.11.09.467903}}

```
