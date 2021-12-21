import numpy as np
from scot2 import *
from evals import *
from utils import *
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import ot 

pca=PCA(n_components=2)

X=np.loadtxt("../data/s1_mapped1.txt")#y=np.genfromtxt("../data/Cheow_methylation.csv", delimiter=",")
y=np.loadtxt("../data/s1_mapped2.txt") #X=np.genfromtxt("../data/Cheow_expression.csv", delimiter=",")

X=unit_normalize(X)
y=unit_normalize(y)

ks=[15]
es=[5e-4]

GW1=[]
GW2=[]
GW3=[]
F=[]
SCOT1=SCOT(X,y)
for k in ks:
	for e in es:
		Xnew, ynew=SCOT1.align(k=k, e=e, balanced=True, verbose=True, normalize=False, XontoY=True)
		gw=SCOT1.coupling
		np.save("simu1_gw_15.npy",gw)
		print(gw.shape)
		print(np.sum(gw))
# 		GW1.append(SCOT1.gwdist)
# 		Cx=get_graph_distance_matrix(Xnew, num_neighbors=k)
# 		Cy=get_graph_distance_matrix(ynew, num_neighbors=k)
# 		gw2=ot.gromov.entropic_gromov_wasserstein2(Cx, Cy, SCOT1.p, SCOT1.q, loss_fun="square_loss", epsilon=e)
# 		GW2.append(gw2)
# 		gw3=ot.gromov.gromov_wasserstein2(Cx, Cy, SCOT1.p, SCOT1.q, loss_fun="square_loss")
# 		GW3.append(gw3)
		print(np.mean(calc_domainAveraged_FOSCTTM(Xnew, ynew)))

# np.save("F.npy", F)
# np.save("GW1.npy", GW1)
# np.save("GW2.npy", GW2)
# np.save("GW3.npy", GW3)

# import matplotlib.pyplot as plt 
# plt.scatter(F, GW1)
# plt.yscale("log")
# plt.show()

# plt.scatter(F, GW2)
# plt.yscale("log")
# plt.show()

# plt.scatter(F, GW3)
# plt.yscale("log")
# plt.show()


# print(GW1)
# print(GW2)
# print(GW3)

