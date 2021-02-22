import numpy as np
from scot2 import *
from evals import *
from utils import *
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

pca=PCA(n_components=2)

X=np.loadtxt("../data/snare_rna.txt")#y=np.genfromtxt("../data/Cheow_methylation.csv", delimiter=",")
y=np.loadtxt("../data/snare_chromatin.txt") #X=np.genfromtxt("../data/Cheow_expression.csv", delimiter=",")
# y=np.loadtxt("../data/snare_rna.txt")
# X=np.loadtxt("../data/snare_chromatin.txt")
# X=unit_normalize(X)
# y=unit_normalize(y)
# x_lab=np.load("../data/sci_rna_labels.npy")
# y_lab=np.load("../data/sci_atac_labels.npy")

types=np.loadtxt("../data/SNAREseq_types.txt")

from collections import Counter
print(Counter(types))
p=[]
p1=[]
q=[]
q1=[]

for i in range(X.shape[0]):
	if types[i]==1:
		p1.append(1.0/(379.0*4.0))
		p.append(379.0/1047.0)
	if types[i]==2:
		p1.append(1.0/(324.0*4.0))
		p.append(324.0/1047.0)
	if types[i]==3:
		p1.append(1.0/(201.0*4.0))
		p.append(201.0/1047.0)
	if types[i]==4:
		p1.append(1.0/(143.0*4.0))
		p.append(143.0/1047.0)

p=p/np.sum(p)
p1=p1/np.sum(p1)

for i in range(y.shape[0]):
	if types[i]==1:
		q1.append(1.0/(379.0*4.0))
		q.append(379.0/1047.0)
	if types[i]==2:
		q1.append(1.0/(324.0*4.0))
		q.append(324.0/1047.0)
	if types[i]==3:
		q1.append(1.0/(201.0*4.0))
		q.append(201.0/1047.0)
	if types[i]==4:
		q1.append(1.0/(143.0*4.0))
		q.append(143.0/1047.0)
q=q/np.sum(q)
q1=q1/np.sum(q1)

print(np.sum(p),np.sum(p1),np.sum(q),np.sum(q))

SCOT1=SCOT(X,y)
SCOT1.p=p1
SCOT1.q=q1
Xnew, ynew=SCOT1.align(k=30,e=5e-4, balanced=True, verbose=True, normalize=True, norm="l2", XontoY=True)
print(np.mean(calc_domainAveraged_FOSCTTM(Xnew,ynew)))
# # print(len(x_lab), len(y_lab))
# # 
# ks=[]
# es=[]
# sils=[]
# for k in [35]:
# 	for e in [5e-3]:
# 		Xnew, ynew= scot(X, y, k, e, mode="connectivity", metric="correlation", XontoY=True, returnCoupling=False, balanced = True)
# 		# Xnew, ynew=SCOT1.align(k=k, e=e, balanced=True, normalize=False) #, norm="l2"
# 		print(np.mean(calc_domainAveraged_FOSCTTM(Xnew,ynew)))
# 		# np.save("SCOTsimu1X.npy",Xnew)
# 		# np.save("SCOTsimu1y.npy",ynew)
# # together=np.concatenate([Xnew, ynew], axis=0)
# # pcad=pca.fit_transform(together)
# # Xnew=pcad[0:177,]
# # ynew=pcad[177:,]


# # plt.scatter(Xnew[:,0], Xnew[:,1], c="black")
# # plt.scatter(ynew[:,0], ynew[:,1], c="red")
# # plt.show()

# np.save("scGEM_alignedExpr.npy", Xnew)
# np.save("scGEM_alignedMethyl.npy", ynew)