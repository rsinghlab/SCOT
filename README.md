
![](assets/SCOT_logo.png)

# Single-Cell alignment using Optimal Transport (SCOT)

SCOT is a python tool for performing unsupervised alignment of single-cell multi-omics datasets. Its methodology is detailed in the pre-print "[Gromov-Wasserstein based optimal transport for aligning single-cell multi-omics data](https://www.biorxiv.org/content/10.1101/2020.04.28.066787v2)"  
For full documentation, please visit https://rsinghlab.github.io/SCOT/ (under construction).  

### SCOT v.0.2.0  
Changes from the previous version:  
**1)** Extensions to cases with differing representation of cell types between domains with unbalanced OT     
We are thankful for [Stephen Zhang](http://zsteve.phatcode.net/about/)'s [contributions](https://github.com/zsteve/SCOT) to our repository with the integration of [unbalanced optimal transport](https://arxiv.org/pdf/2009.04266.pdf) implementation by [Thibault Sejourne](https://github.com/thibsej/unbalanced_gromov_wasserstein)   
**2)** Unsupervised hyperparameter tuning  
When missing correspondence information for performing hyperparameter selection, SCOT uses Gromov-Wasserstein distance as an approximation for alignment quality.  
**3)** Switch to correlation (from Euclidean distance) as a distance metric for kNN graphs  
**4)** Switch to connectivity (from distance) for kNN graph edge annotations  

**Python packages required:**
numpy, sklearn, matplotlib, scipy, cython, POT (note that numpy and cython must be installed prior to POT), torch  

**Folder navivgation:**  
**1) src** contains the source code for SCOT  
**2) data** contains raw data files  
**3) replication** contains jupyter notebooks to replicate results from our paper  
**4) examples** contains short scripts and notebooks to show how to apply SCOT in different scenarios    

**Note:** We are happy to see any work built using or on top of SCOT. However, we ask that you please make sure to give credit in your code if you are using code from this repository.  
Demetci, P. Santorella, R. Sandstede, B., Noble, W. S., Singh, R. 2020. Gromov-Wasserstein based optimal transport for aligning single-cell multi-omics data. bioRxiv.  
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

