---
layout: template
permalink: /tutorial/
---

## Tutorial

#### Installing SCOT

SCOT is developed in Python 3. It depends on a few Python packages, namely: `numpy, cython, scipy, sklearn, matlab, POT` <br>

After you clone our repository locally, you can install these dependencies using the `requirements.txt` file. If you are using `pip`  you can do so by running `pip3 install -r requirements.txt` or `python3 -m pip install -r requirements.txt` on your terminal. If you are using `conda`, you can do so with the command ` conda install --file requirements.txt ` <br>

You can clone the SCOT repository locally in two ways:<br>
**1** If you use `git`, by running `git clone https://github.com/rsinghlab/SCOT.git` on your terminal, or <br>
**2** By navigating to the GitHub repository, clicking on the green `Code` button with the download icon, selecting `Download ZIP` option and then extracting the downloaded compressed folder.  <br>
 
#### Running SCOT

Once you have cloned the SCOT repository and installed the requirements, you will be ready to use it on your own datasets by importing SCOT in a Python script: `from scot import scot`. Note that if your Python script lives elsewhere, you would need to specify the path to scot.py in your local copy of SCOT using `sys`. <br>

SCOT expects datasets to be in `numpy` arrays. 

Once you have read in the datasets, you can initialize the SCOT 

#### Choosing hyperparameters=
If you are not sure 


