---
layout: template
permalink: /tutorial/
---

## Tutorial

#### Installing SCOT
SCOT is developed using Python 3. It depends on a few Python packages, namely: `numpy`, `cython`, `scipy`, `sklearn`, `matlab`, and `POT` <br>

After you clone our repository locally, you can install these dependencies using the `requirements.txt` file.<br>
If you are using `pip`  you can do so by running `pip3 install -r requirements.txt` or `python3 -m pip install -r requirements.txt` on your terminal.<br>
If you are using `conda`, you can do so with the command `conda install --file requirements.txt` <br>

You can clone the SCOT repository locally in one of two ways:<br>
**1** If you use `git`, by running `git clone https://github.com/rsinghlab/SCOT.git` on your terminal, or <br>
**2** By navigating to [our GitHub repository](https://github.com/rsinghlab/SCOT), clicking on the green `Code` button with the download icon, selecting `Download ZIP` option and then extracting the downloaded compressed folder.  <br>
 
#### Running SCOT
Once you have cloned the SCOT repository and installed the requirements, you will be ready to use it on your own datasets by importing SCOT in a Python script:  
`from scot import scot`.  

Note that if your Python script lives elsewhere, you would need to specify the path to scot.py in your local copy of SCOT using `sys`. Example: <br>
```python
import sys
sys.path.insert(1, 'path_to_SCOT')
from scot import SCOT
```
SCOT expects datasets to be in `numpy` arrays. If you have your data in text format, you can read in these using the [`numpy.genfromtxt()`](https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html) or [`numpy.loadtxt()`](https://numpy.org/doc/stable/reference/generated/numpy.loadtxt.html) functions. Example:<br>
```python
import numpy as np 
domain1= np.genfromtxt("path_to_data_file.txt", delimiter="\t") #Change delimiter according to your text file
domain2= np.loadtxt("path_to_data_file2.txt", delimiter="\t") #Same, but with "loadtxt". genfromtxt gives a few more options when loading, e.g. dealing with missing values.
```

If you have `.mtx` data format, which is a common format for single-cell sequencing datasets, you can turn these into `numpy` arrays with the Python package called `scanpy` Example:<br>
```python
import scanpy as sc
my_dataset=sc.read_mtx("datasetFilename.counts.mtx")
my_dataset_npy=my_dataset.X.toarray()
```
Please make sure that the rows in your data matrix/numpy array correspond to samples and columns correspond to genomic features (and transpose your matrix with `numpy.transpose` if needed).<br>

Once you have read in the datasets, you can initialize the SCOT and then run the alignment algorithm on it, which will return the:
```python
import numpy as np
from scot import SCOT

scot_aligner=SCOT(domain1, domain2)
k= 50 # a hyperparameter of the model, determines the number of neighbors to be used in the kNN graph constructed for cells based on sequencing data correlations
e= 1e-3 # another hyperparameter of the model, determines the coefficient of the entropic regularization term
normalize=True #
aligned_domain1, aligned_domain2= scot_aligner.align(k=k, e=epsilon, normalize=normalize)
```
Please take a look at the [examples page](rsinghlab.github.io/SCOT/examples) for Python scripts demonstrating the use of SCOT to align datasets.

#### Choosing hyperparameters
There are two required hyperparameters for performing alignment with SCOT:

| Parameter |       Description     | Default Value | Recommended Range to Try|
| ----------|-----------------------|---------------|-------------------------|
|     k     | Number of neighbors to consider in kNN graphs | 50 | 10 -- n/5, where n is the number of samples (cells) in the smallest dataset |
|     e     | Coefficient of the entropic regularization term in the objective function of OT formulation | 1e-3 | 5e-5 -- 1e-1 |

In general, we have found that the algorithm is fairly robust to the choice of `k` and the parameter `e` makes a larger difference. The larger values of `e` disperses the correspondence probabilities across more samples. If you expect to find 1-to-1 correspondences between samples, err towards smaller values of `e`. 

If you are not sure which hyperparameters to set while running SCOT alignment, you have two options: <br>
**1.** If you have some validation data about the cell-to-cell correspondences between the two domains in your dataset, you can use these for hyperparameter tuning. For this, take a look at [the hyperparameter tuning example script](https://github.com/rsinghlab/SCOT/blob/master/examples/hyperparameterTuning_example.py).<br>
**2.** If you don't have any validation data on correspondences, no worries! You can use [the unsupervised hyperparameter finding procedure](https://rsinghlab.github.io/SCOT/unsupervised/), where we use the Gromov-Wasserstein distance as a proxy for graph distances to check for alignment quality as we sweep through different hyperparameter combinations. <br>

#### Optional parameters for SCOT alignment

align(self, k, e, balanced=True, rho=1e-3, verbose=True, normalize=True, norm="l2", XontoY=True):

| Parameter |       Description     | Default Value | Notes |
| ----------|-----------------------|---------------|-------------------------|-------|
| balanced  | Determines whether to perform balanced or unbalanced optimal tranport | (boolean) `False` | By default, we perform balanced transport. However, if you have a reason to believe there will be severe underrepresentation of at least one cell type in one of the domains in comparison to the other, setting this to False will yield better alignments. Otherwise, keep at True. | 
|    rho    | Coefficient of the Kullback-Leibler relaxation term in unbalanced OT formulation| `5e-2` | Only need this if you are performing unbalanced OT. 
|  verbose  | Determines whether to print transport progress (loss over iterations) while optimizing alignment | (boolean) `True` | 
| normalize | Determines whether to normalize datasets before performing alignment on them. | (boolean) `True` | Empirically, we have found that normalization slightly helps with the resulting alignment quality, so we suggest you keep this `True`.
| normalize | Determines whether to normalize datasets before performing alignment on them. | (boolean) `True` | Empirically, we have found that normalization slightly helps with the resulting alignment quality, so we suggest you keep this `True`.
|    norm   | Defines what sort of normalization will be applied on the datasets. Normalization is always performed column-wise. | `l2` | Empirically, we have found that best alignment results are obtained when the pre-processed (via a dimensionality reduction scheme) real world sequencing datasets are l-2 normalized first. Other options: `zscore`, `l1`, and `max`. 
|  XontoY   | Sets the direction of the barycentric projection. | (boolean) `True` | The direction of the barycentric projection makes very little difference in the quality of the resulting alignment. Generally, it is advisable to project onto the dataset with the more defined clusters. |


#### Attributes of SCOT aligner you can access to:
The `align` function of SCOT returns the two aligned matrices. However, there is additional information you can access:<br> 
- `scot.coupling` will yield the probabilistic correspondence (coupling) matrix. The rows of the matrix will correspond to the samples in the first domain (X) and the columns will correspond to the samples in the second domain (Y). The order of the samples is the same as the order in the input datasets. You can use these correspondence probabilities in your downstream analyses. <br>
- `scot.gwdist` will yield the Gromov-Wasserstein distance between the aligned datasets. This is also internally used when performing fully unsupervised alignment as a proxy for alignment quality. <br>
- `scot.flag` tells whether the optimization procedure for optimal transport has converged. The aligner notifies the user with a printed error message if convergence has failed, but one can also check with this flag. If it has not converged (returns `False`), you might need to set the parameter `e` to a higher value. <br>
- `scot.Xgraph` holds the adjacency matrix of the kNN graph built for the first domain. <br>
- `scot.ygraph` holds the adjacency matrix of the kNN graph built for the second domain. <br>
- `scot.Cx` holds the intra-domain distance matrix for the first domain, computed based on shortest distances on the kNN graph. <br>
- `scot.Cy` holds the intra-domain distance matrix for the second domain, computed based on shortest distances on the kNN graph. <br>
- `scot.p` corresponds to the marginal probability distribution for the samples in the first domain. We use uniform distribution, treating this as a vector of empirical probabilities for each sample. However, if you have some prior information on the marginal probabilities, please change this after initializing SCOT and before running the alignment. <br>
- `scot.q` corresponds to the marginal probability distribution for the samples in the second domain. We use uniform distribution, treating this as a vector of empirical probabilities for each sample. However, if you have some prior information on the marginal probabilities, please change this after initializing SCOT and before running the alignment. 
