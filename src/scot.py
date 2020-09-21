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


def scot(X, y, k, e, mode="distance", metric="minkowski", XontoY=True):
    """
    Given two datasets (X and y) and 
    the hyperparameters (k: number of neighbors to be used in kNN graph construction; and e: eplison value in entropic regularization),
    returns the resulting datasets after transport
    For transport in the opposite direction, set XontoY to False
    """
    # X=ut.zscore_standardize(np.asarray(X))
    # y=ut.zscore_standardize(np.asarray(y))

    Cx=ut.get_graph_distance_matrix(X, k, mode="distance", metric="minkowski") 
    Cy=ut.get_graph_distance_matrix(y, k, mode= "distance", metric="minkowski")
    X_sampleNo= Cx.shape[0]
    y_sampleNo= Cy.shape[0]
    p=ot.unif(X_sampleNo)
    q=ot.unif(y_sampleNo)
    couplingM, log = ot.gromov.entropic_gromov_wasserstein(Cx, Cy, p, q, 'square_loss', epsilon=e, log=True, verbose=True)
    
    # check to make sure GW congerged 
    if (np.isnan(couplingM).any() or np.any(~couplingM.any(axis=1)) or np.any(~couplingM.any(axis=0)) or sum(sum(couplingM)) < .95):
        print("Did not converge. Try increasing the epsilon value. ")
        return None, None
    else:
        if XontoY==True:
            X_transported = ut.transport_data(X,y,couplingM,transposeCoupling=False)
            return X_transported, y
        else:
            y_transported = ut.transport_data(X,y,couplingM,transposeCoupling=True)
            return X, y_transported


    
