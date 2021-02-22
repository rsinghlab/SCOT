from unioncom import UnionCom
import numpy as np
from utils import *
from evals import *
from sklearn.neighbors import KNeighborsClassifier

X=np.loadtxt("scGEM_expression.txt")
y=np.loadtxt("scGEM_methylation.txt")
# X=unit_normalize(X)
# y=unit_normalize(y)
print(X.shape, y.shape)

Xlabel=np.loadtxt("scGEM_typeExpression.txt")
ylabel=np.loadtxt("scGEM_typeMethylation.txt")

uc=UnionCom.UnionCom()
integrated_data=uc.fit_transform(dataset=[X,y])
np.save("scgemUC.npy", integrated_data)
Xn=integrated_data[0]
yn=integrated_data[1]
print(Xn.shape, yn.shape)
print(np.mean(calc_domainAveraged_FOSCTTM(Xn, yn)))
print("LabelTransfer", transfer_accuracy(Xn, yn, Xlabel, ylabel, 5))