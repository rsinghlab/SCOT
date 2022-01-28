
![](assets/SCOT_logo.png)

# Single-Cell alignment using Optimal Transport (SCOT)
* **Note:** We are currently updating this repository with compertmantalized versions of the algorithm due to the new development of v.2.0. Please look for the release of v2.0 version for SCOTv2.

SCOT is a Python tool for performing unsupervised alignment of single-cell multi-omics datasets. Its methodology is detailed in the following two papers:
- SCOT v.1.0: [Gromov-Wasserstein based optimal transport for aligning single-cell multi-omics data](https://www.biorxiv.org/content/10.1101/2020.04.28.066787v2)
- SCOT v.2.0: [Unsupervised integration of single-cell multi-omics datasets with disparities in cell-type representation](https://www.biorxiv.org/content/10.1101/2021.11.09.467903v1)

For full documentation, please visit https://rsinghlab.github.io/SCOT/ (currently being updated).  

## SCOT v.2.0
A few extensions:
1) Alignment with the unbalanced Gromov-Wasserstein optimal transport formulation to handle cell-type representation disparities (Sejourne et al, 2020)
2) Multi-modal alignment by picking the anchor domain based on imputation potential of domain-specific nearest neighbor graphs
3) Different choices for joint embedding/projection

## SCOT v.1.1
A naive extension to multi-modal alignment, where the first dataset in the input as treated as the anchor to align on. 

## SCOT v.1.0
Unsupervised single-cell multi-omic integration with Gromov-Wasserstein optimal transport & a self-tuning heuristic for hyperparameter selection


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
