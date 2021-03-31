"""
Authors: Pinar Demetci, Rebecca Santorella
12 February 2020
Utils for SCOT
"""
import numpy as np
import scipy as sp
from scipy.sparse.csgraph import dijkstra
from scipy.sparse import csr_matrix
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.neighbors import NearestNeighbors, KNeighborsClassifier

def unit_normalize(data, norm="l2", bySample=True):
	"""
	Default norm used is l2-norm. Other options: "l1", and "max"
	If bySample==True, then we independently normalize each sample. If bySample==False, then we independently normalize each feature
	"""
	assert (norm in ["l1","l2","max"]), "Norm argument has to be either one of 'max', 'l1', or 'l2'."

	if bySample==True:
		axis=1
	else:
		axis=0

	return normalize(data, norm=norm, axis=axis) 

def zscore_standardize(data):
	scaler=StandardScaler()
	scaledData=scaler.fit_transform(data)
	return scaledData

def get_spatial_distance_matrix(data, metric="eucledian"):
	Cdata= sp.spatial.distance.cdist(data,data,metric=metric)
	return Cdata/Cdata.max()

def get_graph_distance_matrix(data, num_neighbors, mode="connectivity", metric="correlation"):
	"""
	Compute graph distance matrices on data 
	"""
	assert (mode in ["connectivity", "distance"]), "Norm argument has to be either one of 'connectivity', or 'distance'. "
	if mode=="connectivity":
		include_self=True
	else:
		include_self=False
	graph_data=kneighbors_graph(data, num_neighbors, mode=mode, metric=metric, include_self=include_self)
	shortestPath_data= dijkstra(csgraph= csr_matrix(graph_data), directed=False, return_predecessors=False)
	shortestPath_max= np.nanmax(shortestPath_data[shortestPath_data != np.inf])
	shortestPath_data[shortestPath_data > shortestPath_max] = shortestPath_max
	shortestPath_data=shortestPath_data/shortestPath_data.max()

	return shortestPath_data


def transport_data(source, target, couplingMatrix, transposeCoupling=False):
	"""
	Given: data in the target space, data in the source space, a coupling matrix learned via Gromow-Wasserstein OT
	Returns: source (target) matrix transported onto the target (source)
	"""
	if transposeCoupling == False:
		P = (couplingMatrix.T/couplingMatrix.sum(1)).T
		transported_data= np.matmul(P, target)
	else:
		P = (couplingMatrix/couplingMatrix.sum(0)).T
		transported_data=np.matmul(P, source)
	return transported_data

