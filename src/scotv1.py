"""
Authors: Pinar Demetci, Rebecca Santorella
Principal Investigator: Ritambhara Singh, Ph.D. from Brown University
12 February 2020
Updated: 27 November 2020
SCOT algorithm (version 1): Single Cell alignment using Optimal Transport
Correspondence: pinar_demetci@brown.edu, rebecca_santorella@brown.edu, ritambhara@brown.edu
"""

### Import python packages we depend on:
# For regular matrix operations:
import numpy as np
# For optimal transport operations:
import ot
# For computing graph distances:
from scipy.sparse.csgraph import dijkstra
from scipy.sparse import csr_matrix
from sklearn.neighbors import kneighbors_graph

# For pre-processing, normalization
from sklearn.preprocessing import StandardScaler, normalize


class SCOT(object):
	"""
	SCOT algorithm for unsupervised alignment of single-cell multi-omic data.
	https://www.biorxiv.org/content/10.1101/2020.04.28.066787v2 (original preprint)
	https://www.liebertpub.com/doi/full/10.1089/cmb.2021.0446 (Journal of Computational Biology publication through RECOMB 2021 conference)
	Input: domain1, domain2 in form of numpy arrays/matrices, where the rows correspond to samples and columns correspond to features.
	Returns: aligned domain 1, aligned domain 2 in form of numpy arrays/matrices projected on domain 1
	Example use:
	# Given two numpy matrices, domain1 and domain2, where the rows are cells and columns are different genomic features:
	scot= SCOT(domain1, domain2)
	aligned_domain1, aligned_domain2 = scot.align(k=20, e=1e-3)
	#If you can't pick the parameters k and e, you can try out our unsupervised self-tuning heuristic by running:
	scot= SCOT(domain1, domain2)
	aligned_domain1, aligned_domain2 = scot.align(selfTune=True)
	Required parameters:
	- k: Number of neighbors to be used when constructing kNN graphs. Default= min(min(n_1, n_2), 50), where n_i, for i=1,2 corresponds to the number of samples in the i^th domain.
	- e: Regularization constant for the entropic regularization term in entropic Gromov-Wasserstein optimal transport formulation. Default= 1e-3 
   
	Optional parameters:
	- normalize= Determines whether to normalize input data ahead of alignment. True or False (boolean parameter). Default = True.
	- norm= Determines what sort of normalization to run, "l2", "l1", "max", "zscore". Default="l2" 
	- mode: "connectivity" or "distance". Determines whether to use a connectivity graph (adjacency matrix of 1s/0s based on whether nodes are connected) or a distance graph (adjacency matrix entries weighted by distances between nodes). Default="connectivity"  
	- metric: Sets the metric to use while constructing nearest neighbor graphs. some possible choices are "correlation", "minkowski".  "correlation" is Pearson's correlation and "minkowski" is equivalent to Euclidean distance in its default form (). Default= "correlation". 
	- verbose: Prints loss while optimizing the optimal transport formulation. Default=True
	- XontoY: Determines the direction of barycentric projection. True or False (boolean parameter). If True, projects domain1 onto domain2. If False, projects domain2 onto domain1. Default=True.
	Note: If you want to specify the marginal distributions of the input domains and not use uniform distribution, please set the attributes p and q to the distributions of your choice (for domain 1, and 2, respectively) 
			after initializing a SCOT class instance and before running alignment and set init_marginals=False in .align() parameters
	"""

	def __init__(self, domain1, domain2):

		self.X=domain1
		self.y=domain2

		self.p= None #empirical probability distribution for domain 1 (X)
		self.q= None #empirical probability distribution for domain 2 (y)

		self.Cx=None #intra-domain graph distances for domain 1 (X)
		self.Cy=None #intra-domain graph distances for domain 2 (y)

		self.coupling=None # Coupling matrix that relates domain 1 and domain 2, ..., m
		self.gwdist=None # Gromov-Wasserstein distance between domains after alignment
		self.flag = None # convergence flag

		self.X_aligned=None #aligned datasets to return: domain1
		self.y_aligned=None #aligned datasets to return: domain2

	def init_marginals(self):
		# Without any prior information, we set the probabilities to what we observe empirically: uniform over all observed sample
		self.p= ot.unif(self.X.shape[0])
		self.q = ot.unif(self.y.shape[0])

	def normalize(self, norm="l2", bySample=True):
		assert (norm in ["l1","l2","max", "zscore"]), "Norm argument has to be either one of 'max', 'l1', 'l2' or 'zscore'. If you would like to perform another type of normalization, please give SCOT the normalize data and set the argument normalize=False when running the algorithm."

		if (bySample==True or bySample==None):
			axis=1
		else:
			axis=0

		if norm=="zscore":
			scaler=StandardScaler()
			self.X, self.y=scaler.fit_transform(self.X), scaler.fit_transform(self.y)

		else:
			self.X, self.y =normalize(self.X, norm=norm, axis=axis), normalize(self.y, norm=norm, axis=axis)

	def construct_graph(self, k, mode= "connectivity", metric="correlation"):
		assert (mode in ["connectivity", "distance"]), "Norm argument has to be either one of 'connectivity', or 'distance'. "
		if mode=="connectivity":
			include_self=True
		else:
			include_self=False

		self.Xgraph=kneighbors_graph(self.X, k, mode=mode, metric=metric, include_self=include_self)
		self.ygraph=kneighbors_graph(self.y, k, mode=mode, metric=metric, include_self=include_self)

		return self.Xgraph, self.ygraph

	def init_distances(self):
		# Compute shortest distances
		X_shortestPath=dijkstra(csgraph= csr_matrix(self.Xgraph), directed=False, return_predecessors=False)
		y_shortestPath=dijkstra(csgraph= csr_matrix(self.ygraph), directed=False, return_predecessors=False)

		# Deal with unconnected stuff (infinities):
		X_max=np.nanmax(X_shortestPath[X_shortestPath != np.inf])
		y_max=np.nanmax(y_shortestPath[y_shortestPath != np.inf])
		X_shortestPath[X_shortestPath > X_max] = X_max
		y_shortestPath[y_shortestPath > y_max] = y_max

		# Finnally, normalize the distance matrix:
		self.Cx=X_shortestPath/X_shortestPath.max()
		self.Cy=y_shortestPath/y_shortestPath.max()

		return self.Cx, self.Cy

	def find_correspondences(self, e, verbose=True):
		self.coupling, log= ot.gromov.entropic_gromov_wasserstein(self.Cx, self.Cy, self.p, self.q, loss_fun='square_loss', epsilon=e, log=True, verbose=verbose)
		self.gwdist=log['gw_dist']

		# Check convergence:
		if (np.isnan(self.coupling).any() or np.any(~self.coupling.any(axis=1)) or np.any(~self.coupling.any(axis=0)) or sum(sum(self.coupling)) < .95):
			self.flag=False
		else:
			self.flag=True

		return self.gwdist

	def barycentric_projection(self, XontoY=True):
		if XontoY:
			#Projecting the first domain onto the second domain
			self.y_aligned=self.y
			weights=np.sum(self.coupling, axis = 0)
			self.X_aligned=np.matmul(self.coupling, self.y) / weights[:, None]
		else:
			#Projecting the second domain onto the first domain
			self.X_aligned=self.X
			weights=np.sum(self.coupling, axis = 0)
			self.y_aligned=np.matmul(np.transpose(self.coupling), self.X) / weights[:, None]
		return self.X_aligned, self.y_aligned

	def align(self, k=None, e=1e-3, mode="connectivity", metric="correlation", verbose=True, normalize=True, norm="l2", XontoY=True, selfTune=False, init_marginals=True):
		if normalize:
			self.normalize(norm=norm)
		if init_marginals:
			self.init_marginals()

		if selfTune:
			X_aligned, y_aligned= self.unsupervised_scot()
		else:
			if k==None:
				k=min((int(self.X.shape[0]*0.2), int(self.y.shape[0]*0.2)),50)

			self.construct_graph(k, mode= "connectivity", metric="correlation")
			self.init_distances()
			self.find_correspondences(e=e, verbose=verbose)

			if self.flag==False:
				print("CONVERGENCE ERROR: Optimization procedure runs into numerical errors with the hyperparameters specified. Please try aligning with higher values of epsilon.")
				return
			
			else:
				X_aligned, y_aligned = self.barycentric_projection(XontoY=XontoY)

		self.X_aligned, self.y_aligned=X_aligned, y_aligned
		return self.X_aligned, self.y_aligned

	def search_scot(self, ks, es, all_values = False,  mode= "connectivity", metric="correlation", normalize=True, norm="l2", init_marginals=True):
		'''
		Performs a hyperparameter sweep for given values of k and epsilon
		Default: return the parameters corresponding to the lowest GW distance
		(Optional): return all k, epsilon, and GW values
		'''

		# initialize alignment
		if normalize:
			self.normalize(norm=norm)
		if init_marginals:
			self.init_marginals()

		# Note to self: Incorporate multiprocessing here to speed things up
		# store values of k, epsilon, and gw distance
		total=len(es)*len(ks)
		k_sweep=np.zeros(total)
		e_sweep=np.zeros(total)
		gw_sweep=np.zeros(total)

		gmin = 1
		counter=1

		X_aligned,y_aligned=None, None
		e_best,k_best=None, None
		# search in k first to reduce graph computation
		for k in ks:
			self.construct_graph(k, mode= mode, metric=metric)
			self.init_distances()
			for e in es:
				print(counter, "/", total)
				print("Aligning k: ",k, " and e: ",e)
				# run alignment / optimize correspondence matrix:
				self.find_correspondences(e=e, verbose=False)
				# save values
				if self.flag:
					if all_values:
						k_sweep[counter]=k
						e_sweep[counter]=e
						gw_sweep[counter] = self.gwdist

						print(self.gwdist)
					# save the alignment if it is lower
					if self.gwdist < gmin:
						X_aligned, y_aligned = self.barycentric_projection()
						gmin =self.gwdist
						e_best, k_best= e, k
					counter = counter + 1
		   
		if all_values:
			# return alignment and all values
			return X_aligned, y_aligned, gw_sweep, k_sweep, e_sweep
		else:
			# return  alignment and the parameters corresponding to the lowest GW distance
			return X_aligned, y_aligned, gmin, k_best, e_best


	def unsupervised_scot(self, normalize=False, norm='l2'):
		'''
		Unsupervised hyperparameter tuning algorithm to find an alignment
		by using the GW distance as a measure of alignment
		'''

		# use k = 20% of # sample or k = 50 if dataset is large
		n = min(self.X.shape[0], self.y.shape[0])
		k_start = min(n // 5, 50)

		num_eps = 12
		num_k = 5

		# define search space
		es = np.logspace(-1, -3, num_eps)
		if ( n > 250):
			ks = np.linspace(20, 100, num_k)
		else:
			ks = np.linspace(n//20, n//6, num_k)
		ks = ks.astype(int)
		
		# search parameter space
		X_aligned, y_aligned, g_best, k_best, e_best = self.search_scot(ks, es, all_values=False, normalize=normalize, norm=norm, init_marginals=False)

		print("Alignment completed. Hyperparameters selected from the unsupervised hyperparameter sweep are: %d for number of neighbors k and %f for epsilon" %(k_best, e_best))

		return X_aligned, 