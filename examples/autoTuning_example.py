"""
@Author: Pinar Demetci, Rebecca Santorella 2020
SCOT hyperparameter tuning example script
"""

import sys
# We change the working directory because this is where the source code we want to import exists.
# This is not needed if you put the source code (scot.py) in the same directory as your scripts. 
sys.path.insert(1, '../src/')
from scot import *
import evals

### Read and normalize the data:
X=np.load("../data/scatac_feat.npy")
y=np.load("../data/scrna_feat.npy")

# initialize SCOT object
scot=SCOT(X, y)

# perform unsupervised hyperparameter search
X_aligned, y_aligned = scot.align(selfTune=True)
# Hyperparameter combinations picked automatically prints after the unsupervised sweep 

FOSCTTM_unsupervised=np.mean(evals.calc_domainAveraged_FOSCTTM(X_aligned, y_aligned))
print("The average FOSCTTM measure for the unsupervised alignment is: ", FOSCTTM_unsupervised)