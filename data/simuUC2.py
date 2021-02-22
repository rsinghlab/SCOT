from unioncom import UnionCom
import numpy as np
from utils import *
from evals import *

X=np.loadtxt("s2_mapped1.txt")
y=np.loadtxt("s2_mapped2.txt")
# X=zscore_standardize(X)
# y=zscore_standardize(y)
# print(X.shape, y.shape)

Xlabel=np.loadtxt("s2_type1.txt")
ylabel=np.loadtxt("s2_type2.txt")

uc=UnionCom.UnionCom()
integrated_data=uc.fit_transform(dataset=[X,y])
np.save("simu2UC_norm.npy", integrated_data)
Xn=integrated_data[0]
yn=integrated_data[1]
print(Xn.shape, yn.shape)
print(np.mean(calc_domainAveraged_FOSCTTM(Xn, yn)))
print("LabelTransfer", transfer_accuracy(Xn, yn, Xlabel, ylabel, 5))