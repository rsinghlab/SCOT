"""
@Author: Pinar Demetci, Rebecca Santorella 2020
SCOT hyperparameter tuning example script
"""

import os
import sys
sys.path.insert(1, '../src/')
import utils as ut
import evals as evals
import scot2 as sc
import numpy as np

### Change working directory to /data in order to import the data
os.chdir("../data/")

### Read and normalize the data:
X=np.load("scatac_feat.npy")
y=np.load("scrna_feat.npy")

# initialize SCOT object
scot=sc.SCOT(X, y)

# perform unsupervised hyperparameter search
X_aligned, y_aligned, k_best, e_best = scot.unsupervised_scot(XontoY=False)
FOSCTTM_y=np.mean(evals.calc_domainAveraged_FOSCTTM(X_aligned, y_aligned))

print("Best performing setting in the y onto X projection direction is:")
print("k= ",k_best, " epsilon= ",e_best)
print("with an average FOSCTTM measure of: ", FOSCTTM_y)

