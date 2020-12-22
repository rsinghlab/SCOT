---
layout: template
permalink: /tutorial/
---

## Tutorial

#### Installing SCOT

SCOT is developed in Python 3. It depends on a few Python packages, namely: `numpy`, `cython`, `scipy`, `sklearn`, `matlab`, and `POT` <br>

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

SCOT expects datasets to be in `numpy` arrays. If you have your data in text format, you can read in these using the [`numpy.genfromtxt()`]() or [`numpy.loadtxt()`]() functions. Example:<br>
```python
import numpy as np 
domain1= np.genfromtxt("path_to_data_file.txt", delimiter="\t") #Change delimiter according to your text file
domain2= np.loadtxt("path_to_data_file2.txt", delimiter="\t") #Same, but with "loadtxt". genfromtxt gives a few more options when loading, e.g. dealing with missing values.
```

If you have `.mtx` data format, which is a common format for single-cell sequencing datasets, you can turn these into `numpy` arrays with the Python package called `..`

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
If you are not sure which hyperparameters to set while running SCOT alignment, you have two options: <br>
**1.** If you have some validation data about the cell-to-cell correspondences between the two domains, you can use these for hyperparameter tuning. For this, take a look at [the hyperparameter tuning example script]().
**2.** If you have 

