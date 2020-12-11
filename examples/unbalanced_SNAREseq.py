"""
@Author: Pinar Demetci, Rebecca Santorella 2020
SCOT hyperparameter tuning example script
"""

import sys
sys.path.insert(1, '../src/')
import utils as ut
import evals as evals
import scot2 as sc
import numpy as np

X=np.load("../data/scatac_feat.npy")
y=np.load("../data/scrna_feat.npy")

scot=sc.SCOT(X, y)
X_aligned, y_aligned = scot.align(k=50, e=0.0005, balanced=True, rho=5e-2, verbose=True, normalize=True, norm="l2", XontoY=True)

print(scot.gwdist)
print(scot.coupling)
print(np.sum(scot.coupling))

print(X_aligned.shape, y_aligned.shape)
print(np.mean(evals.calc_domainAveraged_FOSCTTM(X_aligned, y_aligned))) 
