"""
Author: Pinar Demetci
Principal Investigator: Ritambhara Singh, Ph.D. from Brown University
08 August 2021
Updated: 23 February 2023
SCOTv2 algorithm: Single Cell alignment using Optimal Transport version 2
Correspondence: pinar_demetci@brown.edu, ritambhara@brown.edu
"""

### Import python packages we depend on:
import numpy as np
import torch
import ot
import scipy
# For computing graph distances:
from scipy.sparse.csgraph import dijkstra
from scipy.sparse import csr_matrix
from sklearn.neighbors import kneighbors_graph

# For pre-processing, normalization
from sklearn.preprocessing import StandardScaler, normalize


class SCOTv2(object):
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

	def __init__(self, data):

		assert type(data)==list and len(data)>=2, "As input, SCOTv2 requires a list, containing at least two numpy arrays to be aligned.  \
				Each numpy array/matrix corresponds to a dataset, with samples (cells) in rows and features (latent representations or genomic features) in columns. \
				We recommend using latent representations (e.g. principal components for RNA-seq and topics - via cisTopic- for ATAC-seq/Methyl-seq)."
		self.data=data
		self.marginals=[] # Holds the empirical probability distributions over samples in each dataset
		self.graphs=[] # Holds graphs per dataset
		self.graphDists=[] # Holds intra-domain graph distances for each input dataset
		self.couplings=[] # Holds coupling matrices
		self.gwdists=[] # Gromov-Wasserstein distances between domains after alignment
		self.flags = [] # Holds alignment convergence flags (booleans: True/False)

		self.aligned_data=[]

	def _init_marginals(self):
		# Without any prior information, we set the probabilities to what we observe empirically: uniform over all observed sample
		for i in range(len(self.data)):
			num_cells=self.data[i].shape[0]
			marginalDist=torch.ones(num_cells)/num_cells
			self.marginals.append(marginalDist)
		return self.marginals

	def _normalize(self, norm="l2", bySample=True):
		assert (norm in ["l1","l2","max", "zscore"]), "Norm argument has to be either one of 'max', 'l1', 'l2' or 'zscore'.\
		 If you would like to perform another type of normalization, please give SCOT the normalized data and set the argument 'normalize=False' when running the algorithm. \
		 We have found l2 normalization to empirically perform better with single-cell sequencing datasets, including when using latent representations. "

		for i in range(len(self.data)):
			if norm=="zscore":
				scaler=StandardScaler()
				self.data[i]=scaler.fit_transform(self.data[i])
			else:
				if (bySample==True or bySample==None):
					axis=1
				else:
					axis=0
				self.data[i] =normalize(self.data[i], norm=norm, axis=axis)
		return self.data # Normalized data

	def construct_graph(self, k=20, mode= "connectivity", metric="correlation"):
		assert (mode in ["connectivity", "distance"]), "Norm argument has to be either one of 'connectivity', or 'distance'. "
		if mode=="connectivity":
			include_self=True
		else:
			include_self=False

		for i in range(len(self.data)):
			self.graphs.append(kneighbors_graph(self.data[i], n_neighbors=k, mode=mode, metric=metric, include_self=include_self))

		return self.graphs

	def init_graph_distances(self):
		for i in range(len(self.data)):
			# Compute shortest distances
			shortestPath=dijkstra(csgraph= csr_matrix(self.graphs[i]), directed=False, return_predecessors=False)
			# Deal with unconnected stuff (infinities):
			Max_dist=np.nanmax(shortestPath[shortestPath != np.inf])
			shortestPath[shortestPath > Max_dist] = Max_dist
			# Finnally, normalize the distance matrix:
			self.graphDists.append(shortestPath/shortestPath.max())

		return self.graphDists

	def _exp_sinkhorn_solver(self, ecost, u, v,a,b, mass, eps, rho, rho2, nits_sinkhorn, tol_sinkhorn):
			"""
			Parameters
			----------
			- ecost: torch.Tensor of size [size_X, size_Y]
					 Exponential kernel generated from the local cost based on the current coupling.  
			- u: torch.Tensor of size [size_X[0]].
				 First dual potential defined on X.
			- v: torch.Tensor of size [size_Y[0]].
				 Second dual potential defined on Y. 
			- mass: torch.Tensor of size [1]. 
					Mass of the current coupling.
			- nits_sinkhorn: int. 
							 Maximum number of iterations to update Sinkhorn potentials in inner loop.
			- tol_sinkhorn: float
							Tolerance on convergence of Sinkhorn potentials.

			Returns
			----------
			u: torch.Tensor of size [size_X[0]]
			   First dual potential of Sinkhorn algorithm
			v: torch.Tensor of size [size_Y[0]]
			   Second dual potential of Sinkhorn algorithm
			logpi: torch.Tensor of size [size_X, size_Y]
				   Optimal transport plan in log-space.
			"""
			# Initialize potentials by finding best translation
			if u is None or v is None:
				u, v = torch.ones_like(a), torch.ones_like(b)
			k = (a * u ** (-eps / rho)).sum()+ (b * v ** (-eps / rho)).sum()
			k = k / (2 * (u[:, None] * v[None, :] * ecost *a[:, None] * b[None, :]).sum())
			z = (0.5 * mass * eps) / (2.0 + 0.5 * (eps / rho) + 0.5 * (eps / rho2))
			k = k ** z
			u,v= u * k, v * k

			# perform Sinkhorn updates in LSE form
			for j in range(nits_sinkhorn):
				u_prev = u.clone()
				v = torch.einsum("ij,i->j", ecost, a * u) ** (-1.0 / (1.0 + eps / rho))
				u = torch.einsum("ij,j->i", ecost, b * v) ** (-1.0 / (1.0 + eps / rho2))
				if (u.log() - u_prev.log()).abs().max().item() * eps < tol_sinkhorn:
					break
			pi = u[:, None] * v[None, :] * ecost * a[:, None] * b[None, :]
			return u, v, pi

	def exp_unbalanced_gw(self,a, dx, b, dy, eps=0.01, rho=1.0, rho2=None, nits_plan=3000, tol_plan=1e-6, nits_sinkhorn=3000, tol_sinkhorn=1e-6):
		if rho2 is None:
			rho2 = rho #KL divergence coefficient doesn't have to be the same for both couplings. 
					   #But, to keep #hyperparameters low, we default to using the same coefficient. 
					   #Someone else playing with our code could assign a rho2 different than rho, though.

		# Initialize the coupling and local costs
		pi= a[:, None]* b[None, :] / (a.sum() * b.sum()).sqrt()
		pi_prev = torch.zeros_like(pi)
		up, vp = None, None

		for i in range(nits_plan):
			pi_prev = pi.clone()
			mp = pi.sum()

			#Compute the current local cost:
			distxy = torch.einsum("ij,kj->ik", dx, torch.einsum("kl,jl->kj", dy, pi))
			kl_pi = torch.sum(pi * (pi / (a[:, None] * b[None, :]) + 1e-10).log())
			mu, nu = torch.sum(pi, dim=1), torch.sum(pi, dim=0)
			distxx = torch.einsum("ij,j->i", dx ** 2, mu)
			distyy = torch.einsum("kl,l->k", dy ** 2, nu)
			lcost = (distxx[:, None] + distyy[None, :] - 2 * distxy) + eps * kl_pi
			if rho < float("Inf"):
				lcost = (lcost+ rho* torch.sum(mu * (mu / a + 1e-10).log()))
			if rho2 < float("Inf"):
				lcost = (lcost+ rho2* torch.sum(nu * (nu / b + 1e-10).log()))
			ecost = (-lcost /(mp * eps)).exp()

			if (i%10)==0:
				print("Unbalanced GW step:", i)
			#compute the coupling via sinkhorn
			up, vp, pi = self._exp_sinkhorn_solver(ecost, up, vp, a, b, mp, eps, rho, rho2,nits_sinkhorn, tol_sinkhorn)
			
			flag=True
			if torch.any(torch.isnan(pi)):
				flag=False

			pi = (mp / pi.sum()).sqrt() * pi
			if (pi - pi_prev).abs().max().item() < tol_plan:
				break
		return pi, flag

	def find_correspondences(self, normalize=True, norm="l2", bySample=True, k=20, mode= "connectivity", metric="correlation",  eps=0.01, rho=1.0, rho2=None):
		# Normalize 
		if normalize:
			self._normalize(norm=norm, bySample=bySample)
		# Initialize inputs for (unbalanced) Gromov-Wasserstein optimal transport:
		self._init_marginals()
		print("computing intra-domain graph distances")
		self.construct_graph(k=k, mode=mode, metric=metric)
		self.init_graph_distances()
		# Run pairwise dataset alignments:
		for i in range(len(self.data)-1):
			print("running pairwise dataset alignments")
			a,b =torch.Tensor(self.marginals[0]), torch.Tensor(self.marginals[i+1])
			dx, dy= torch.Tensor(self.graphDists[0]), torch.Tensor(self.graphDists[i+1])
			coupling, flag=self.exp_unbalanced_gw(a, dx, b, dy, eps=eps, rho=rho, rho2=rho2, nits_plan=3000, tol_plan=1e-6, nits_sinkhorn=3000, tol_sinkhorn=1e-6)
			self.couplings.append(coupling)
			self.flags.append(flag)
			if flag==False:
					raise Exception(
					f"Solver got NaN plan with params (eps, rho, rho2) "
					f" = {eps, rho, rho2}. Try increasing argument eps")
		return self.couplings

	def barycentric_projection(self):
		aligned_datasets=[self.data[0]]
		for i in range(0,len(self.couplings)):
			coupling=np.transpose(self.couplings[i].numpy())
			weights=np.sum(coupling, axis = 1)
			projected_data=np.matmul((coupling/ weights[:, None]), self.data[0])
			aligned_datasets.append(projected_data)
		return aligned_datasets

	def coembed_datasets(self, Lambda=1.0, out_dim=10):
		"""
		Co-embeds datasets in a shared space.
		Implementation is based on Cao et al 2022 (Pamona)
		"""
		n_datasets = len(self.data)
		H0 = []
		L = []
		for i in range(n_datasets-1):
			self.couplings[i] = self.couplings[i]*np.shape(self.data[i])[0]

		for i in range(n_datasets):    
			graph_data = self.graphs[i] + self.graphs[i].T.multiply(self.graphs[i].T > self.graphs[i]) - \
				self.graphs[i].multiply(self.graphs[i].T > self.graphs[i])
			W = np.array(graph_data.todense())
			index_pos = np.where(W>0)
			W[index_pos] = 1/W[index_pos] 
			D = np.diag(np.dot(W, np.ones(np.shape(W)[1])))
			L.append(D - W)

		Sigma_x = []
		Sigma_y = []
		for i in range(n_datasets-1):
			Sigma_y.append(np.diag(np.dot(np.transpose(np.ones(np.shape(self.coupling[i])[0])), self.coupling[i])))
			Sigma_x.append(np.diag(np.dot(self.coupling[i], np.ones(np.shape(self.coupling[i])[1]))))

		S_xy = coupling[0]
		S_xx = L[0] + Lambda*Sigma_x[0]
		S_yy = L[-1] +Lambda*Sigma_y[0]
		for i in range(1, n_datasets-1):
			S_xy = np.vstack((S_xy, self.coupling[i]))
			S_xx = block_diag(S_xx, L[i] + Lambda*Sigma_x[i])
			S_yy = S_yy + Lambda*Sigma_y[i]

		v, Q = np.linalg.eig(S_xx)
		v = v + 1e-12   
		V = np.diag(v**(-0.5))
		H_x = np.dot(Q, np.dot(V, np.transpose(Q)))

		v, Q = np.linalg.eig(S_yy)
		v = v + 1e-12      
		V = np.diag(v**(-0.5))
		H_y = np.dot(Q, np.dot(V, np.transpose(Q)))

		H = np.dot(H_x, np.dot(S_xy, H_y))
		U, sigma, V = np.linalg.svd(H)

		num = [0]
		for i in range(n_datasets-1):
			num.append(num[i]+len(data[i]))

		U, V = U[:,:output_dim], np.transpose(V)[:,:output_dim]

		fx = np.dot(H_x, U)
		fy = np.dot(H_y, V)

		integrated_data = []
		for i in range(n_datasets-1):
			integrated_data.append(fx[num[i]:num[i+1]])

		integrated_data.append(fy)

		return integrated_data

	def align(self,normalize=True, norm="l2", bySample=True, k=20, mode= "connectivity", metric="correlation",  eps=0.01, rho=1.0, rho2=None, projMethod="embedding", Lambda=1.0, out_dim=10):
		assert projMethod in ["embedding", "barycentric"], "The input to the parameter 'projMethod' needs to be one of \
								'embedding' (if co-embedding them in a new shared space) or 'barycentric' (if using barycentric projection)"
		self.find_correspondences(normalize=normalize, norm=norm, bySample=bySample, k=k, mode=mode, metric=metric,  eps=eps, rho=rho, rho2=rho2)
		print("FLAGS", self.flags)
		if projMethod=="embedding":
			integrated_data=self.coembed_datasets(Lambda=Lambda, out_dim=out_dim)
		else:
			integrated_data=self.barycentric_projection()
		self.integrated_data=integrated_data
		return integrated_data
	
# X=np.load("../data/SNARE/SNAREseq_atac_feat.npy")[0:1000,:]
# Y=np.load("../data/SNARE/SNAREseq_rna_feat.npy")
# print(X.shape, Y.shape)
# SCOT=SCOTv2([Y,X])
# aligned_datasets=SCOT.align(normalize=True, k=50, eps=0.005, rho=0.1, projMethod="barycentric")
# print(len(aligned_datasets))
# print(aligned_datasets[0].shape)
# print(aligned_datasets[1].shape)
# # np.save("aligned_atac.npy", aligned_datasets[1])
# np.save("aligned_rna.npy", aligned_datasets[0])


