"""
Authors: Pinar Demetci, Rebecca Santorella
Principal Investigator: Ritambhara Singh, Ph.D. from Brown University
12 February 2020
Updated: 27 November 2020
SCOT algorithm: Single Cell alignment using Optimal Transport
Correspondence: pinar_demetci@brown.edu, rebecca_santorella@brown.edu, ritambhara@brown.edu
"""

### Import python packages we depend on:
# For regular matrix operations:
import numpy as np
# For optimal transport operations:
import ot
from ot.unbalanced import sinkhorn_unbalanced
from ot.gromov import init_matrix, gwloss, gwggrad
# For computing graph distances:
from sklearn.neighbors import kneighbors_graph
from scipy.sparse.csgraph import dijkstra
# For pre-processing, normalization
from sklearn.preprocessing import StandardScaler, normalize
# For convergence errors and parameter warnings:
import sys
import warnings


class SCOT(object):
    """
    SCOT algorithm for unsupervised alignment of single-cell multi-omic data.
    https://www.biorxiv.org/content/10.1101/2020.04.28.066787v2
    
    Example use:
    scot= SCOT(domain1, domain2)
    aligned_domain1, aligned_domain2= scot.align(k, e, balanced=True, rho=None, verbose=True, normalize=True, norm="l2", XontoY=True)

    Input: domain1, domain2, in form of numpy arrays/matrices, where the rows correspond to samples and columns correspond to features.
    Returns: aligned domain 1, aligned domain 2, in form of numpy arrays/matrices.

    Parameters: 
    k: Number of neighbors to be used when constructing kNN graphs. Default is min(50, 0.2n), where n is the number of samples in the smallest domain
    e: Regularization constant for the entropic regularization term in entropic Gromov-Wasserstein optimal transport formulation. Default: 1e-4. We recommend users to search a grid between 5e-4 to 1e-2
    balanced: If you believe there will be a significant underrepresentation/overrepresentation of certain cell types in one of the domains you attempt to align, set this to False. When set to False, it performs unbalanced optimal transport to account for underrepresentation. Default=True. 
    rho: Only needs to be set if using unbalanced OT (if balanced is set to False). Defines the regularization constant for KL relaxation term in unbalanced optimal transport. Default = 5e-2. Ideal value defines on the extent of underrepresentation of cell types between the domains (more unbalanced might want more relaxation)
    verbose: Prints loss when optimizing the optimal transport formulation. Default=True
    normalize: When set to True, normalizes the input domains before performing alingment, otherwise skips normalization. Default= True
    norm: Describes what type of data normalization to use. Available options: "l2", "l1", "max", "zscore". Default= "l2". We have found l2 normalization yields better empirical results with real world single-cell sequencing data.
    XontoY: Describes the direction of barcentric projection used for alignment. When set to True, projects domain 1 onto domain 2. False does opposite. Direction of projection makes little difference in alignment quality. Default= True.
    """
    def __init__(self, domain1, domain2, normalize=True):
        self.X=domain1
        self.y=domain2

        self.p= None #empirical probability distribution for domain 1 (X)
        self.q= None #empirical probability distribution for domain 2 (y)

        self.Xgraph=None #kNN graph of domain 1 (X)
        self.ygraphh=None #kNN graph of domain 2 (y)
        self.Cx=None #intra-domain graph distances for domain 1 (X)
        self.Cy=None #intra-domain graph distances for domain 2 (y)

        self.coupling=None # Coupling matrix that relates domain 1 and domain 2. Entries describes the probability of correspondence between the samples in domain 1 (rows) and domain 2 (columns)
        self.gwdist=None # Gromov-Wasserstein distance between domains after alignment. Can be used as a proxy for alignment quality
        self.flag = None # convergence flag

    def init_marginals(self):
        self.p= ot.unif(self.X.shape[0]) # Without any prior information, we set the probabilities to what we observe empirically: uniform over all observed samples
        self.q= ot.unif(self.y.shape[0]) # Without any prior information, we set the probabilities to what we observe empirically: uniform over all observed samples

    def normalize(self, norm="l2", bySample=True):
        assert (norm in ["l1","l2","max", "zscore"]), "Norm argument has to be either one of 'max', 'l1', 'l2' or 'zscore'. If you would like to perform another type of normalization, please give SCOT the normalize data and set the argument normalize=False when running the algorithm."

        if (bySample==True or bySample==None):
            axis=1
        else:
            axis=0

        if norm=="zscore":
            scaler=StandardScaler()
            self.X=scaler.fit_transform(self.X)
            self.y=scaler.fit_transform(self.y)
        else:
            self.X=normalize(self.X, norm=norm, axis=axis)
            self.y=normalize(self.y, norm=norm, axis=axis)

    def build_kNN(self, k, mode="connectivity", metric="correlation"):
        """
        Helper function: Builds kNN graphs for each domain
        To be used in intradomain distance matrix computations
        """
        
        # Check that inputs are legal
        assert (mode in ["connectivity", "distance"]), "Mode argument has to be either one of 'connectivity' or 'distance'. "
        assert(k <= min(self.X.shape[0], self.y.shape[0])), "Please set the argument k (for the number of neighbors in kNN graphs) to something no larger than the number of samples in the domain with the fewest samples"

        if mode=="connectivity":
            include_self=True
        else:
            include_self=False
        
        # Build graph for the two domains
        self.Xgraph= kneighbors_graph(self.X, k, mode=mode, metric=metric, include_self=include_self)
        self.ygraph= kneighbors_graph(self.y, k, mode=mode, metric=metric, include_self=include_self)

    def compute_graphDistances(self):
        # Compute shortest paths on the kNN graphs for intra-domain distances:
        self.Cx=dijkstra(csgraph= self.Xgraph, directed=False, return_predecessors=False)
        self.Cy=dijkstra(csgraph=  self.ygraph, directed=False, return_predecessors=False)

        # Checking for illegal values: if infinite distances exist (e.g. not a connected graph), correct them to maximum finite distance on the graph:
        X_maxShortest= np.nanmax(self.Cx[self.Cx != np.inf])
        self.Cx[self.Cx > X_maxShortest] = X_maxShortest

        y_maxShortest= np.nanmax(self.Cy[self.Cy != np.inf])
        self.Cy[self.Cy > y_maxShortest] = y_maxShortest

        # Normalize intra-domain distances based on maximum distance:
        self.Cx=np.asarray(self.Cx/self.Cx.max(), dtype=np.float64)
        self.Cy=np.asarray(self.Cy/self.Cy.max(), dtype=np.float64)

    def unbalanced_entropic_gromov_wasserstein(self, e, rho, loss_fun="square_loss", max_iter=1000, tol=1e-6, verbose=True):
        """
        Helper function: unbalanced Gromov-Wasserstein OT for when there is cell type unbalance between domains
        Used when balanced=False in align() function.
        Adapted from POT package using ot.unbalanced.sinkhorn() and ot.gromov.entropic_gromov_wasserstein()

        Parameters:
        e: Regularization constant for the entropic regularization term in entropic Gromov-Wasserstein optimal transport formulation.
        rho: Regularization constant for KL relaxation term in unbalanced Sinkhorn-Knopp iterations. 

        Returns: Coupling matrix and log dictionary where error and GW distance from optimization have been logged
        """
        coupling = np.outer(self.p, self.q)  # Initialize the coupling matrix
        constC, hCx, hCy = init_matrix(self.Cx, self.Cy, self.p, self.q, loss_fun)

        cpt = 0
        err = 1    

        log = {'err': []}

        while (err > tol and cpt < max_iter):
            couplingPrev = coupling
            # compute the gradient
            grad = gwggrad(constC, hCx, hCy, coupling)
            coupling = sinkhorn_unbalanced(self.p, self.q, grad, e, rho, method='sinkhorn', numItermax=max_iter, stopThr=tol, verbose=verbose, log=False)
            
            if cpt % 10 == 0:
                # we can speed up the process by checking for the error only all the 10th iterations
                err = np.linalg.norm(coupling - couplingPrev)
                log['err'].append(err)

                if verbose:
                    if cpt % 200 == 0:
                        print('{:5s}|{:12s}'.format(
                            'It.', 'Err') + '\n' + '-' * 19)
                print('{:5d}|{:8e}|'.format(cpt, err))
            cpt += 1

        log['gw_dist'] = gwloss(constC, hCx, hCy, coupling)
        return coupling, log

    def find_correspondences(self, e, balanced=True, rho=5e-2, verbose=True):
        self.flag = True
        if balanced:
            self.coupling, log = ot.gromov.entropic_gromov_wasserstein(self.Cx, self.Cy, self.p, self.q, 'square_loss', epsilon=e, log=True, verbose=verbose)
            self.gwdist= log['gw_dist']

            if (np.isnan(self.coupling).any() or np.isinf(self.coupling).any() or np.sum(self.coupling) < .98):
                #sys.exit("Alignment algorithm did not converge. This is very likely due to low e values (the epsilon parameter set for entropic regularization constant). Please try again with a higher e value. We recommend not going below 5e-4.")
                self.flag = False

        else:
            self.coupling, log = self.unbalanced_entropic_gromov_wasserstein(e=e, rho= rho, loss_fun="square_loss", max_iter=1000, tol=1e-6, verbose=True)
            self.gwdist= log['gw_dist']

            if (np.isnan(self.coupling).any() or np.isinf(self.coupling).any()): # Note: It's possible for unbalanced OT to converge and return a coupling matrix where sum does not add up to one. So we cannot check for sum of coupling elements for coupling in this case. 
                #sys.exit("Alignment algorithm did not converge. This is very likely due to low e values (the epsilon parameter set for entropic regularization constant). Please try again with a higher e value. We recommend not going below 5e-4.")
                self.flag = False

    
    def barycentric_projection(self, XontoY=True):

        if XontoY:
            y_aligned=self.y
            X_weights = np.sum(self.coupling, axis = 1)
            X_aligned=np.matmul(self.coupling, self.y) / X_weights[:,None]
        else:
            X_aligned=self.X
            y_weights = np.sum(self.coupling, axis = 0)
            y_aligned=np.matmul(np.transpose(self.coupling), self.X) / y_weights[:,None]
        return X_aligned, y_aligned

    def align(self, k, e, balanced=True, rho=1e-3, verbose=True, normalize=True, norm="l2", XontoY=True):
        if normalize:
            self.normalize(norm=norm)
        self.init_marginals()
        self.build_kNN(k)
        self.compute_graphDistances()
        self.find_correspondences(e=e, balanced=balanced,  rho=rho, verbose=verbose)
        X_aligned, y_aligned = self.barycentric_projection(XontoY)
        return X_aligned, y_aligned

    def search_scot(self, ks, es, all_values = False): 
        '''
        Performs a hyperparameter sweep for given values of k and epsilon
        Default: return the parameters corresponding to the lowest GW distance
        (Optional): return all k, epsilon, and GW values
        '''

        self.init_marginals()

        # store values of k, epsilon, and gw distance 
        k_plot=[]
        e_plot=[]
        g_plot=[]

        total=len(es)*len(ks)
        counter=0
        
        # search in k first to reduce graph computation
        for k in ks:
            self.build_kNN(k)
            self.compute_graphDistances()

            for e in es:

                counter+=1
                if (counter % 10 == 0):
                    print(str(counter)+"/"+str(total))

                # run scot
                self.find_correspondences(e=e, balanced=True)

                if self.flag:   
                    g_plot.append(self.gwdist)
                    k_plot.append(k)
                    e_plot.append(e)          
       
        if all_values:
                # return all values     
            return g_plot, k_plot, e_plot
        else:
            # return the parameters corresponding to the lowest gromov-wasserstein distance
            gmin=np.amin(g_plot)
            gminI=np.argmin(g_plot)
            e_best = e_plot[gminI]
            k_best = k_plot[gminI]
            return k_best, e_best

    def unsupervised_scot(self, normalize=True, norm='l2', XontoY=True):
        '''
        Unsupervised hyperparameter tuning algorithm to find an alignment 
        by using the GW distance as a measure of alignment
        '''

        if normalize:
            self.normalize(norm=norm)

        # use k = 20% of # sample or k = 50 if dataset is large 
        n = min(self.X.shape[0], self.y.shape[0]) 
        k_best = min(n // 5, 50)

        # first fix k and find the best epsilon (6 runs)
        es = np.logspace(-2, -3, 6)
        g1, k1, e1 = self.search_scot([k_best], es, all_values = True)

        # save the best epsilon from that search
        gmin = np.min(g1)
        gminI=np.argmin(g1)
        e_best = e1[gminI]

        # fix that epsilon and vary k (4 runs)
        if ( n > 250):
            ks = np.linspace(20, 100, 4)
        else:
            ks = np.linspace(n//20, n//6, 4)
        ks = ks.astype(int)
        g2, k2, e2 = self.search_scot(ks, [e_best], all_values = True)

        # save the best k from that search 
        gminI=np.argmin(g2)
        if (g2[gminI] < gmin):
            gmin = g2[gminI]
            k_best = k2[gminI]

        # now use that k and epsilon to do a more refined grid search (10 runs)
        scale = np.log10(e_best)
        eps_refined = np.logspace(scale + .25, scale - .25, 5)

        ks_refined = np.linspace( max(5, k_best - 5), min(n//2, k_best + 5), 2)    
        ks_refined = ks_refined.astype(int)
        g3, k3, e3 = self.search_scot(ks_refined, eps_refined, all_values = True)

        # find the best parameter set from all runs
        gminI=np.argmin(g3)
        if (g3[gminI] < gmin):
            gmin = g3[gminI]
            k_best = k3[gminI]
            e_best = e3[gminI]

        # get alignment with these parameters
        self.build_kNN(k_best)
        self.compute_graphDistances()
        self.find_correspondences(e=e_best, balanced=True)
        X_aligned, y_aligned = self.barycentric_projection(XontoY)
        
        return X_aligned, y_aligned, k_best, e_best

