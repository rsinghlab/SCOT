"""
Author: Pinar Demetci
29 September 2020
SCOT on GPU
"""
import numpy as np
import ot
from ot.utils import dist, UndefinedParameter
from ot.optim import cg
from ot.gromov import init_matrix, gwggrad, gwloss
import src.utils as ut
from ot.gpu import bregman

def entropic_gromov_wasserstein(Cx, Cy, p, q, loss_fun, epsilon,  max_iter=1000, tol=1e-9, verbose=False, log=False):
	"""
	Carries out Gromov-Wasserstein OT by running Sinkhorn-Knopp iterations on the GPU
	"""

    C1 = np.asarray(C1, dtype=np.float64)
    C2 = np.asarray(C2, dtype=np.float64)

    T = np.outer(p, q)  # Initialization

    constC, hC1, hC2 = init_matrix(C1, C2, p, q, loss_fun)

    cpt = 0
    err = 1

    if log:
        log = {'err': []}

    while (err > tol and cpt < max_iter):

        Tprev = T

        # compute the gradient
        tens = gwggrad(constC, hC1, hC2, T)

        T = bregman.sinkhorn_knopp(p, q, tens, epsilon)

        if cpt % 10 == 0:
            # we can speed up the process by checking for the error only all
            # the 10th iterations
            err = np.linalg.norm(T - Tprev)

            if log:
                log['err'].append(err)

            if verbose:
                if cpt % 200 == 0:
                    print('{:5s}|{:12s}'.format(
                        'It.', 'Err') + '\n' + '-' * 19)
                print('{:5d}|{:8e}|'.format(cpt, err))

        cpt += 1

    if log:
        log['gw_dist'] = gwloss(constC, hC1, hC2, T)
        return T, log
    else:
        return T


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
	couplingM, log = ot.gromov.entropic_gromov_wasserstein(Cx, Cy, p, q, loss='square_loss', epsilon=e, log=True, verbose=True)

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
