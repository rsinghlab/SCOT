"""
@Author: Pinar Demetci, 2020
SCOT hyperparameter tuning script
"""

import os
import numpy as np
import src.utils as ut
import src.evals as evals
from src.scot import *

### Change working directory to /data in order to import the data
os.chdir("data/")

### Read and normalize the data:
X=np.load("scatac_feat.npy")
y=np.load("scrna_feat.npy")
X=ut.unit_normalize(X)
y=ut.unit_normalize(y)

### Set the grid of hyperparameters to try:
es=[5e-4, 7e-4, 1e-3, 3e-3,  5e-3, 7e-3, 1e-2, 3e-2, 5e-2, 7e-2, 1e-1, 3e-1, 5e-1, 7e-1]
ks=[5,10,20,30,40,50,60,70, 80,90,100]
total=len(es)*len(ks) #total number of hyperparameters

### Initialize lists for recording:
all_ks=[]
all_es=[]
all_FOSCTTM_X=[]
all_FOSCTTM_y=[]
all_GWdist=[]



counter=1
for e in es:
	for k in ks:
		print(counter," of ", total)
		couplingM, log=scot(X, y, k, e, mode="connectivity", metric="correlation", returnCoupling=True)	

		# Check convergence: 
		converged=True
		if (np.isnan(couplingM).any() or np.any(~couplingM.any(axis=1)) or np.any(~couplingM.any(axis=0)) or sum(sum(couplingM)) < .95):
			converged=False

		# Get metrics:
		gwDist=log['gw_dist']

		# Project X onto Y and calculate mean FOSCTTM:
		X_transported = ut.transport_data(X,y,couplingM,transposeCoupling=False)
		FOSCTTM_X=np.mean(evals.calc_domainAveraged_FOSCTTM(X_transported,y))

		# Project X onto Y and calculate mean FOSCTTM:
		y_transported = ut.transport_data(X,y,couplingM,transposeCoupling=True)
		FOSCTTM_y=np.mean(evals.calc_domainAveraged_FOSCTTM(X,y_transported))

		# Record parameters and metrics:
		all_ks.append(k)
		all_es.append(e)
		if converged==True:
			all_FOSCTTM_X.append(FOSCTTM_X)
			all_FOSCTTM_y.append(FOSCTTM_y)
			all_GWdist.append(gwDist)

		else: # So that convergence issues don't mislead us to a numerically unstable hyperparameter combination:
			all_FOSCTTM_X.append(1)
			all_FOSCTTM_y.append(1)
			all_GWdist.append(10)
			
		counter+=1

# Choose the best performing hyperparameter setting in both directions:
index_X = np.where(all_FOSCTTM_X==np.amin(all_FOSCTTM_X))[0][0]
print("index_X", index_X)
print("Best performing setting in the X onto y projection direction is:")
print("k= ",all_ks[index_X], " epsilon= ",all_es[index_X])
print("with an average FOSCTTM measure of: ", all_FOSCTTM_X[index_X])

index_y = np.where(all_FOSCTTM_y==np.amin(all_FOSCTTM_y))[0][0]
print("Best performing setting in the y onto X projection direction is:")
print("k= ",all_ks[index_y], " epsilon= ",all_es[index_y])
print("with an average FOSCTTM measure of: ", all_FOSCTTM_y[index_y])

# Save these lists (if desired) so we can investigate correlations, trends, spread of metrics etc:
# (could also turn into a matrix or dictionary and save that way so there aren't too many files)
np.save("snare_tuning_es.npy", all_es)
np.save("snare_tuning_ks.npy", all_ks)
np.save("snare_tuning_FOSCTTM_X.npy", all_FOSCTTM_X)
np.save("snare_tuning_FOSCTTM_y.npy", all_FOSCTTM_y)
np.save("snare_tuning_GWdist.npy", all_GWdist)
