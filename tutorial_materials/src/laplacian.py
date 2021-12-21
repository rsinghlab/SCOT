import numpy as np 

D1=np.asarray([[2,0,0,0,0,0],[0,3,0,0,0,0],[0,0,3,0,0,0],[0,0,0,2,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
A1=np.asarray([[0,1,1,0,0,0],[1,0,1,0,0,1],[1,1,0,1,0,0],[0,0,1,0,1,0],[0,0,0,1,0,0],[0,1,0,0,0,0]])

D2=np.asarray([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,2,0,0,0],[0,0,0,3,0,0],[0,0,0,0,3,0],[0,0,0,0,0,1]])
A2=np.asarray([[0,1,0,0,0,0],[1,0,0,1,0,0],[0,0,0,1,1,0],[0,1,1,0,1,0],[0,0,1,1,0,1],[0,0,0,0,1,0]])

L1=D1-A1
L2=D2-A2

# print(L1)
w1,v1= np.linalg.eig(L1)
# print(w1)
# print(v1)

# print(L2)
w2,v2= np.linalg.eig(L2)
# print(w2)
# print(v2)

print(w1*v1)
print(w2*v2)
# u1,s1,vh1= np.linalg.svd(L1)
# u2,s2,vh2= np.linalg.svd(L2)

# print(u1)
# print(u2)
# print(s1)
# print(s2)
# print(vh1)
# print(vh2)