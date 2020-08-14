"""
Author: Pinar Demetci, Rebecca Santorella
12 February 2020
Utils for SCOT
"""
import ot
from ot.bregman import sinkhorn
from ot.utils import dist, UndefinedParameter
from ot.optim import cg
from ot.gromov import init_matrix, gwggrad, gwloss
import evaluation_metrics as em
import utils as ut

def stabilized_entropic_gromov_wasserstein(C1, C2, p, q, loss_fun, epsilon,
                                max_iter = 1000, tol=1e-9, verbose=False, log=False):


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

        T = sinkhorn(p, q, tens, epsilon, method = 'sinkhorn_stabilized')

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

def scot(X, y, k, e, XontoY=True):
    """
    Given two datasets (X and y) and 
    the hyperparameters (k: number of neighbors to be used in kNN graph construction; and e: eplison value in entropic regularization),
    returns the resulting datasets after transport
    For transport in the opposite direction, set XontoY to False
    """
    X=ut.zscore_standardize(np.asarray(X))
    y=ut.zscore_standardize(np.asarray(y))

    Cx=ut.get_graph_distance_matrix(X, k, mode="distance") 
    Cy=ut.get_graph_distance_matrix(y, k, mode= "distance")
    X_sampleNo= Cx.shape[0]
    y_sampleNo= Cy.shape[0]
    p=ot.unif(X_sampleNo)
    q=ot.unif(y_sampleNo)
    couplingM, log = stabilized_entropic_gromov_wasserstein(Cx, Cy, p, q, 'square_loss', epsilon=e, log=True, verbose=True)

    if XontoY==True:
        X_transported = ut.transport_data(X,y,couplingM,transposeCoupling=False)
        return X_transported, y
    else:
        y_transported = ut.transport_data(X,y,couplingM,transposeCoupling=True)
        return X, y_transported


    
