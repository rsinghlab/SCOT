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


def scot(X, y, k, e, mode="distance", metric="minkowski", XontoY=True, returnCoupling=False):
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
	Cy=ut.get_graph_distance_matrix(y, k, mode=mode, metric=metric)

	# Initialize marginal distributions over data:
	X_sampleNo= Cx.shape[0]
	y_sampleNo= Cy.shape[0]
	p=ot.unif(X_sampleNo)
	q=ot.unif(y_sampleNo)

	# Perform optimization to get the coupling matrix between domains:
	couplingM, log = ot.gromov.entropic_gromov_wasserstein(Cx, Cy, p, q, 'square_loss', epsilon=e, log=True, verbose=True)

	# check to make sure GW congerged, if not, warn the user with an error statement
	converged=True #initialize the convergence flag
	if (np.isnan(couplingM).any() or np.any(~couplingM.any(axis=1)) or np.any(~couplingM.any(axis=0)) or sum(sum(couplingM)) < .95):
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

    if plot:
        return g_plot, k_plot, e_plot
    else:
        return k_best, e_best

# find the best alignment by gromov-wasserstein distance
def unsupervised_scot(X,y, XontoY=True):
    
    # first fix k=50 and find the best epsilon (20 runs)
    es = np.logspace(-1, -4, 20)
    k, eps = sweep_scot(X,y,[50], es)

    # now use that epsilon to find the best k (10 runs)
    ks = [10,20,30,40,50,60,70,80,90,100]
    k, eps = sweep_scot(X,y,ks, [eps])

    # now use that k and epsilon to do a more refined grid search (30 runs)
    scale = np.log10(eps)
    eps_refined = np.logspace(scale + .25, scale - .25, 6)
    if k > 10:
        ks_refined = [k - 10, k - 5, k, k + 5, k + 10]
    else:
        ks_refined = [k - 5, k, k + 5, k + 10, k + 15]
    k_best, e_best = search_scot(X, y, ks_refined, eps_refined)

    # run soct with these parameters
    X_t, y_t = scot(X, y, k_best, e_best, XontoY)
    
    return X_t, y_t, k_best, e_best


