import numpy as np
from scot import *
from evals import *
from utils import *
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

pca=PCA(n_components=2)

y=np.genfromtxt("../data/Cheow_methylation.csv", delimiter=",")
X=np.genfromtxt("../data/Cheow_expression.csv", delimiter=",")
# y=np.loadtxt("../data/snare_rna.txt")
# X=np.loadtxt("../data/snare_chromatin.txt")
X=unit_normalize(X)
y=unit_normalize(y)
# x_lab=np.load("../data/sci_rna_labels.npy")
# y_lab=np.load("../data/sci_atac_labels.npy")

# print(len(x_lab), len(y_lab))
# SCOT1=SCOT(X,y)
ks=[]
es=[]
sils=[]
for k in [35]:
	for e in [5e-3]:
		Xnew, ynew= scot(X, y, k, e, mode="connectivity", metric="correlation", XontoY=True, returnCoupling=False, balanced = True)
		# Xnew, ynew=SCOT1.align(k=k, e=e, balanced=True, normalize=False) #, norm="l2"
		print(np.mean(calc_domainAveraged_FOSCTTM(Xnew,ynew)))
		# np.save("SCOTsimu1X.npy",Xnew)
		# np.save("SCOTsimu1y.npy",ynew)
# together=np.concatenate([Xnew, ynew], axis=0)
# pcad=pca.fit_transform(together)
# Xnew=pcad[0:177,]
# ynew=pcad[177:,]


# plt.scatter(Xnew[:,0], Xnew[:,1], c="black")
# plt.scatter(ynew[:,0], ynew[:,1], c="red")
# plt.show()

np.save("scGEM_alignedExpr.npy", Xnew)
np.save("scGEM_alignedMethyl.npy", ynew)