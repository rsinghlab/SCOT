"""
Author: Pinar Demetci, Rebecca Santorella
12 February 2020
Utils for SCOT
"""
import numpy as np
import ot
from ot.bregman import sinkhorn
from ot.utils import dist, UndefinedParameter
from ot.optim import cg
from ot.gromov import init_matrix, gwggrad, gwloss
import src.utils as ut
import sys

# We use the unbalanced Gromov-Wasserstein solver (utilising PyTorch) written by Thibault Sejourne 
# https://github.com/thibsej/unbalanced_gromov_wasserstein
# Clone the above repo and update below to the corresponding path
sys.path.insert(0, "/home/zsteve/analysis/unbalanced_gromov_wasserstein/solver")
from tlb_kl_sinkhorn_solver import TLBSinkhornSolver
import torch

def scot(X, y, k, e, rho = 1, mode="connectivity", metric="correlation", XontoY=True, returnCoupling=False, balanced = True):
	"""
	Given two datasets (X and y) and 
	the hyperparameters (k: number of neighbors to be used in kNN graph construction; and e: eplison value in entropic regularization),
	returns the resulting datasets after transport
	For transport in the opposite direction, set XontoY to False
	"""
	## I think we should let users choose what type of normalization they want to use since that might matter in comparisons to other methods
	## e.g. MMD-MA uses l-2 norm sometimes and they might want to use that for fair comparison. So commented out the part below -Pinar
	# X=ut.zscore_standardize(np.asarray(X))
	# y=ut.zscore_standardize(np.asarray(y))

	# Construct the kNN graphs
	Cx=ut.get_graph_distance_matrix(X, k, mode=mode, metric=metric) 
	Cy=ut.get_graph_distance_matrix(y, k, mode=mode, metric=metric);

	# Initialize marginal distributions over data:
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
	if (np.isnan(couplingM).any() or np.any(~couplingM.any(axis=1)) or np.any(~couplingM.any(axis=0))): # or sum(sum(couplingM)) < .95):
		print("Did not converge. Try increasing the epsilon value. ")
		converged=False

	# If the user wants to get the coupling matrix and the optimization log at the end of this and investigate it or perform some projection themselves,
	# allow them (useful for hyperparameter tuning with projections in both directions) :
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
# runs scot for given values of k and epsilon
# by default returns the parameters corresponding to the lowest gromov-wasserstein distance
# (optional) returns all data points for plotting
def search_scot(X,y, ks, es, plot_values = False): 

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
            gw, log = ot.gromov.entropic_gromov_wasserstein(Cx, Cy, p, q, 'square_loss', epsilon = e, log=True, verbose=False, max_iter = 200)

            if (np.isnan(gw).any() or np.any(~gw.any(axis=1)) or np.any(~gw.any(axis=0)) or sum(sum(gw)) < .95):
                print("Did not converge")
            else:   
                g_plot.append(log["gw_dist"])
                k_plot.append(k)
                e_plot.append(e)          

    # find the parameters corresponding to the lowest gromov-wasserstein distance
    gmin=np.amin(g_plot)
    gminI=np.argmin(g_plot)
    e_best = e_plot[gminI]
    k_best = k_plot[gminI]
    print("Best result with GW distance is when e and k are:", e_best, k_best, " with lowest GW dist:", gmin)    

    if plot_values:
        return g_plot, k_plot, e_plot
    else:
        return k_best, e_best

# find the best alignment by gromov-wasserstein distance
def unsupervised_scot(X,y, XontoY=True):

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

    print("Lowest GW distance is ", gmin, " for epsilon = ", e_best, " and k = ", k_best)

    # run soct with these parameters
    X_t, y_t = scot(X, y, k_best, e_best, XontoY = XontoY)
    
    return X_t, y_t, k_best, e_best
