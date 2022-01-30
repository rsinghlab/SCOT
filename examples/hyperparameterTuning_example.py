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
X=np.genfromtxt("../data/scGEM_expression.csv", delimiter=",")
y=np.genfromtxt("../data/scGEM_expression.csv", delimiter=",")

### Set the grid of hyperparameters to try:
es=[5e-4, 7e-4, 1e-3, 3e-3,  5e-3, 7e-3, 1e-2, 3e-2, 5e-2, 7e-2, 1e-1, 3e-1, 5e-1, 7e-1]
ks=[5, 10,20,30,40,50,60,70, 80,90,100]
total=len(es)*len(ks) #total number of hyperparameters

### Initialize lists for recording:
all_ks=[]
all_es=[]
all_FOSCTTM_X=[]
all_FOSCTTM_y=[]
all_GWdist=[]

# initialize SCOT object
scot=SCOT(X, y)

counter=1
for e in es:
    for k in ks:
        print(counter," of ", total)
        X_aligned, y_aligned = scot.align(k, e, normalize = True, norm="l2", XontoY=True)    
        
        # check convergence
        if scot.flag:
            # Record parameters and metrics:
            all_ks.append(k)
            all_es.append(e)
            
            # Use X onto Y projection to calculate mean FOSCTTM:
            FOSCTTM_X=np.mean(evals.calc_domainAveraged_FOSCTTM(X_aligned, y_aligned))

            # Use Y onto Y projection to calculate mean FOSCTTM:
            X_aligned2, y_aligned2 = scot.barycentric_projection(XontoY=False)
            FOSCTTM_y=np.mean(evals.calc_domainAveraged_FOSCTTM(X_aligned2, y_aligned2))

            all_FOSCTTM_X.append(FOSCTTM_X)
            all_FOSCTTM_y.append(FOSCTTM_y)
            all_GWdist.append(scot.gwdist)

        else: # If scot.flag is False, it means we had convergence issues for these hyperparameter combination
            # We do the following to record large FOSCTTM and GWdist values so that convergence issues don't mislead us to a numerically unstable hyperparameter combination
            # Alternatively, you can comment out everything here and just write "pass" as to not record any numerically unstable configurations at all:

            #pass

            all_FOSCTTM_X.append(1) #Largest FOSCTTM possible
            all_FOSCTTM_y.append(1) #Largest FOSCTTM possible
            all_GWdist.append(10) #A very large GW distance
            # Record parameters and metrics:
            all_ks.append(k)
            all_es.append(e)

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
np.save("../data/scGEM_tuning_es.npy", all_es)
np.save("../data/scGEM_tuning_ks.npy", all_ks)
np.save("../data/scGEM_tuning_FOSCTTM_X.npy", all_FOSCTTM_X)
np.save("../data/scGEM_tuning_FOSCTTM_y.npy", all_FOSCTTM_y)
np.save("../data/scGEM_tuning_GWdist.npy", all_GWdist)
