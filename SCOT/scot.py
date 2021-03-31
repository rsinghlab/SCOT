"""
Author: Pinar Demetci, Rebecca Santorella
12 February 2020
Updated: 27 November 2020
SCOT algorithm: Single Cell alignment using Optimal Transport
"""
import numpy as np
import ot
from ot.bregman import sinkhorn
from ot.utils import dist, UndefinedParameter
from ot.optim import cg
from ot.gromov import init_matrix, gwggrad, gwloss
import .utils as ut
import sys

from unbalancedgw import TLBSinkhornSolver
import torch

def scot(X, y, k, e, rho = 1, mode="connectivity", metric="correlation", XontoY=True, returnCoupling=False, balanced = True):
	"""
	Given two datasets (X and y) and 
	the hyperparameters (k: number of neighbors to be used in kNN graph construction; and e: eplison value in entropic regularization),
	returns the resulting datasets after transport
	For transport in the opposite direction, set XontoY to False
	"""
	
	# Construct the kNN graphs
	Cx=ut.get_graph_distance_matrix(X, k, mode=mode, metric=metric) 
	Cy=ut.get_graph_distance_matrix(y, k, mode=mode, metric=metric)

	# Initialize uniform marginal distributions over data:
	X_sampleNo= Cx.shape[0]
	y_sampleNo= Cy.shape[0]
	p=ot.unif(X_sampleNo)
	q=ot.unif(y_sampleNo)

	# Perform optimization to get the coupling matrix between domains:
	if balanced:
		couplingM, log = ot.gromov.entropic_gromov_wasserstein(Cx, Cy, p, q, 'square_loss', epsilon=e, log=True, verbose=True)
	else:
		solver = TLBSinkhornSolver(nits=1000, nits_sinkhorn=2500, gradient=False, tol=1e-3, tol_sinkhorn=1e-3)
		couplingM, _ = solver.tlb_sinkhorn(torch.Tensor(p).cuda(), torch.Tensor(Cx).cuda(), torch.Tensor(q).cuda(), torch.Tensor(Cy).cuda(), rho=rho*0.5*(Cx.mean() + Cy.mean()), eps=e, init=None)
		couplingM = couplingM.cpu().numpy()
		log = None

	# check to make sure GW congerged, if not, warn the user with an error statement
	converged=True #initialize the convergence flag
	if (np.isnan(couplingM).any() or np.any(~couplingM.any(axis=1)) or np.any(~couplingM.any(axis=0))):
		print("Did not converge. Try increasing the epsilon value. ")
		converged=False

	# Allow the user to return the coupling if desired
	if returnCoupling==True:
		return couplingM, log

	# Otherwise perform barycentric projection in the desired direction and return the aligned matrices
	else:
		if converged==False: #except if the convergence failed, just return None, None.
			return None, None
		if XontoY==True:
			X_transported = ut.transport_data(X,y,couplingM,transposeCoupling=False)
			return X_transported, y
		else:
			y_transported = ut.transport_data(X,y,couplingM,transposeCoupling=True)
			return X, y_transported
		
def search_scot(X,y, ks, es, plot_values = False): 
    '''
    Performs a hyperparameter sweep for given values of k and epsilon
    Default: return the parameters corresponding to the lowest GW distance
    (Optional): return all k, epsilon, and GW values
    '''
    X_sampleNo= X.shape[0]
    y_sampleNo= y.shape[0]
    p=ot.unif(X_sampleNo)
    q=ot.unif(y_sampleNo)

    # store values of k, epsilon, and gw distance 
    k_plot=[]
    e_plot=[]
    g_plot=[]

    total=len(es)*len(ks)
    counter=0
    
    # search in k first to reduce graph computation
    for k in ks:
        Cx=ut.get_graph_distance_matrix(X, k, mode="connectivity", metric="correlation") 
        Cy=ut.get_graph_distance_matrix(y, k, mode="connectivity", metric="correlation")

        for e in es:

            counter+=1
            if (counter % 10 == 0):
                print(str(counter)+"/"+str(total))

            # run scot
            gw, log = ot.gromov.entropic_gromov_wasserstein(Cx, Cy, p, q, 'square_loss', epsilon = e, log=True, verbose=False)

            if (np.isnan(gw).any() or np.any(~gw.any(axis=1)) or np.any(~gw.any(axis=0))):
                print("Did not converge")
            else:   
                g_plot.append(log["gw_dist"])
                k_plot.append(k)
                e_plot.append(e)          
   
    if plot_values:
	# return all values	
        return g_plot, k_plot, e_plot
    else:
        # return the parameters corresponding to the lowest gromov-wasserstein distance
        gmin=np.amin(g_plot)
        gminI=np.argmin(g_plot)
        e_best = e_plot[gminI]
        k_best = k_plot[gminI]
        return k_best, e_best

def unsupervised_scot(X,y, XontoY=True):
    '''
    Unsupervised hyperparameter tuning algorithm to find an alignment 
    by using the GW distance as a measure of alignment
    '''
    # use k = 20% of # sample or k = 50 if dataset is large 
    n = min(X.shape[0], y.shape[0]) 
    k_best = min(n // 5, 50)

    # first fix k and find the best epsilon (6 runs)
    es = np.logspace(-2, -3, 6)
    g1, k1, e1 = search_scot(X,y,[k_best], es, plot_values = True)

    # save the best epsilon from that search
    gmin = np.min(g1)
    gminI=np.argmin(g1)
    e_best = e1[gminI]

    # fix that epsilon and vary k (4 runs)
    if ( n > 250):
        ks = np.linspace(20, 100, 4)
    else:
        ks = np.linspace(X.shape[0]//20, X.shape[0]//6, 4)
    ks = ks.astype(int)
    g2, k2, e2 = search_scot(X,y,ks, [e_best], plot_values = True)

    # save the best k from that search 
    gminI=np.argmin(g2)
    if (g2[gminI] < gmin):
        gmin = g2[gminI]
        k_best = k2[gminI]

    # now use that k and epsilon to do a more refined grid search (10 runs)
    scale = np.log10(e_best)
    eps_refined = np.logspace(scale + .25, scale - .25, 5)

    ks_refined = np.linspace( max(5, k_best - 5), min(X.shape[0]//2, k_best + 5), 2)    
    ks_refined = ks_refined.astype(int)
    g3, k3, e3 = search_scot(X, y, ks_refined, eps_refined, plot_values = True)

    # find the best parameter set from all runs
    gminI=np.argmin(g3)
    if (g3[gminI] < gmin):
        gmin = g3[gminI]
        k_best = k3[gminI]
        e_best = e3[gminI]

    # run soct with these parameters
    X_t, y_t = scot(X, y, k_best, e_best, XontoY = XontoY)
    
    return X_t, y_t, k_best, e_best
