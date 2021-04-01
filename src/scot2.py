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
import utils as ut


class SCOT(object):
    """
    SCOT algorithm for unsupervised alignment of single-cell multi-omic data.
    https://www.biorxiv.org/content/10.1101/2020.04.28.066787v2

    Example use:
    scot= SCOT(domain1, [domain2, domain3, ... domaink])
    aligned_domain1, aligned_domain2tok = scot.align(k, e, balanced=True, rho=None, verbose=True, normalize=True, norm="l2", XontoY=True)
    Input: domain1, domain2,... domain k in form of numpy arrays/matrices, where the rows correspond to samples and columns correspond to features.
    Returns: aligned domain 1, [aligned domain 2, ... aligned domain k] in form of numpy arrays/matrices projected on domain 1
    Parameters:
    k: Number of neighbors to be used when constructing kNN graphs. D
    e: Regularization constant for the entropic regularization term in entropic Gromov-Wasserstein optimal transport formulation. 
    balanced: If you believe there will be a significant underrepresentation/overrepresentation of certain cell types in one of the domains you attempt to align, set this to False. When set to False, it performs unbalanced optimal transport to account for underrepresentation. Default=True.
    rho: Only needs to be set if using unbalanced OT (if balanced is set to False). Defines the regularization constant for KL relaxation term in unbalanced optimal transport. Default = 5e-2. Ideal value defines on the extent of underrepresentation of cell types between the domains (more unbalanced might want more relaxation)
    verbose: Prints loss when optimizing the optimal transport formulation. Default=True
    normalize: When set to True, normalizes the input domains before performing alingment, otherwise skips normalization. Default= True
    norm: Describes what type of data normalization to use. Available options: "l2", "l1", "max", "zscore". Default= "l2". We have found l2 normalization yields better empirical results with real world single-cell sequencing data.
    """
    def __init__(self, domain1, domain2, normalize=True):
        self.X=domain1
        if type(domain2) == list:
            self.y=domain2
        else:
            self.y = [domain2]
        self.num_y = len(self.y)

        self.p= None #empirical probability distribution for domain 1 (X)
        self.q= None #empirical probability distribution for domain 2 (y)

        self.Cx=None #intra-domain graph distances for domain 1 (X)
        self.Cy=None #intra-domain graph distances for domain 2 (y)

        self.coupling=None # Coupling matrix that relates domain 1 and domain 2, ..., m
        self.gwdist=None # Gromov-Wasserstein distance between domains after alignment
        self.flag = None # convergence flag

    # Without any prior information, we set the probabilities to what we observe empirically: uniform over all observed sample
    def init_marginals(self):
        self.p= ot.unif(self.X.shape[0])
        self.q = [ot.unif(self.y[iy].shape[0]) for iy in range(self.num_y)]

    def init_coupling(self):
        self.coupling = [np.outer(self.p, self.q[iy]) for iy in range(self.num_y)]

    def normalize(self, norm="l2", bySample=True):
        assert (norm in ["l1","l2","max", "zscore"]), "Norm argument has to be either one of 'max', 'l1', 'l2' or 'zscore'. If you would like to perform another type of normalization, please give SCOT the normalize data and set the argument normalize=False when running the algorithm."

        if (bySample==True or bySample==None):
            axis=1
        else:
            axis=0

        if norm=="zscore":
            scaler=StandardScaler()
            self.X=scaler.fit_transform(self.X)
            self.y=[scaler.fit_transform(self.y[iy]) for iy in range(self.num_y)]
        else:
            self.X=normalize(self.X, norm=norm, axis=axis)
            self.y = [normalize(self.y[iy], norm=norm, axis=axis) for iy in range(self.num_y)]

    def init_distances(self, k, mode="connectivity", metric="correlation"):
        self.Cx = []
        self.Cy = []
        if type(k) != list:
            k = [k for i in range(self.num_y)]
        self.Cx = [ut.get_graph_distance_matrix(self.X, k[iy],  mode, metric) for iy in range(self.num_y)]
        self.Cy = [ut.get_graph_distance_matrix(self.y[iy], k[iy], mode, metric) for iy in range(self.num_y)]

    def entropic_gromov_wasserstein(self, iy, loss_fun, epsilon, rho=0.05, max_iter=1000, tol=1e-9, verbose=False, balanced = True):
        """
        Adapted from POT package using ot.unbalanced.sinkhorn() and ot.gromov.entropic_gromov_wasserstein()

        Returns the gromov-wasserstein transport

        If balanced = True, the function solves the following optimization problem:

        .. math::
            GW = arg\min_T \sum_{i,j,k,l} L(Cx_{i,k},Cy_{j,l})*T_{i,j}*T_{k,l}-\epsilon(H(T))

            s.t. T 1 = p

                 T^T 1= q

                 T\geq 0

        If balanced = False, the function solves the following optimization problem:

        .. math::
            GW = arg\min_T \sum_{i,j,k,l} L(Cx_{i,k},Cy_{j,l})*T_{i,j}*T_{k,l}-\epsilon(H(T))
                                + \rho KL(T 1 || p) + \rho KL(T^T 1 || q)

            s.t.  T\geq 0

        Where :
        - Cx : Metric cost matrix in the source space
        - Cy : Metric cost matrix in the target space
        - p  : distribution in the source space
        - q  : distribution in the target space
        - L  : loss function to account for the misfit between the similarity matrices
        - H  : entropy

        Parameters
        ----------
        Cx : ndarray, shape (ns, ns)
            Metric cost matrix in the source space
        Cy : num_y length list of ndarray, shape (nt, nt)
            Metric costfr matrix in the target space
        p :  ndarray, shape (ns,)
            Distribution in the source space
        q : num_y length list of ndarray, shape (nt,)
            Distribution in the target space
        loss_fun :  string
            Loss function used for the solver either 'square_loss' or 'kl_loss'
        epsilon : float
            Regularization term >0
        rho : float
            Regularization term > 0
        max_iter : int, optional
            Max number of iterations
        tol : float, optional
            Stop threshold on error (>0)
        verbose : bool, optional
            Print information along iterations
        balanced: bool, optional
            Compute balanced or unbalanced GW

        Returns
        -------
        T : ndarray, shape (ns, nt)
            Optimal coupling between the two spaces

        References
        ----------
        .. [12] Peyré, Gabriel, Marco Cuturi, and Justin Solomon,
            "Gromov-Wasserstein averaging of kernel and distance matrices."
            International Conference on Machine Learning (ICML). 2016.

        """
        self.Cx[iy] = np.asarray(self.Cx[iy], dtype=np.float64)
        self.Cy[iy] =  np.asarray(self.Cy[iy], dtype=np.float64)
        T = self.coupling[iy]
        constC, hCx, hCy = ot.gromov.init_matrix(self.Cx[iy], self.Cy[iy], self.p, self.q[iy], loss_fun)

        cpt = 0
        err = 1

        log = {'err': []}

        while (err > tol and cpt < max_iter):
            Tprev = T

            # compute the gradient
            tens =ot.gromov.gwggrad(constC, hCx, hCy, T)
            if balanced:
                T = ot.bregman.sinkhorn(self.p, self.q[iy], tens, epsilon)
            else:
                T = ot.unbalanced.sinkhorn_unbalanced(self.p, self.q[iy], tens, epsilon,
                                rho, method='sinkhorn', numItermax=max_iter, stopThr=tol, verbose=verbose, log=False)

            if cpt % 10 == 0:
                # we can speed up the process by checking for the error only all
                # the 10th iterations
                err = np.linalg.norm(T - Tprev)
                log['err'] = err

                if verbose:
                    if cpt % 200 == 0:
                        print('{:5s}|{:12s}'.format(
                            'It.', 'Err') + '\n' + '-' * 19)
                    print('{:5d}|{:8e}|'.format(cpt, err))

            cpt += 1

        log['num_iter'] = cpt
        log['gw_dist'] = ot.gromov.gwloss(constC, hCx, hCy, T) 
        return T, log

    def find_correspondences(self, e, balanced=True, rho=5e-2, verbose=True):
        self.flag = [True for iy in range(self.num_y)]
        self.gwdist = [False for iy in range(self.num_y)]
        if type(e) != list:
            e = [e for i in range(self.num_y)]
        if type(rho) != list:
            rho = [rho for i in range(self.num_y)]
        for i in range(self.num_y):
            m, log = self.entropic_gromov_wasserstein(i, loss_fun="square_loss", epsilon=e[i], rho=rho[i],
                                verbose=verbose, balanced = balanced)
            # check convergence
            if ( np.isnan(m).any() or np.isinf(m.any()) ):
                self.flag[i] = False
            else:
                self.coupling[i] = m
                self.gwdist[i] = log['gw_dist']
        return 

    def barycentric_projection(self):
        X_aligned=self.X
        y_aligned = []
        for i in range(self.num_y):
            weights = np.sum(self.coupling[i], axis = 0)
            y_aligned.append( np.matmul(np.transpose(self.coupling[i]), self.X) / weights[:, None])
        return X_aligned, y_aligned

    def align(self, k, e, balanced=True, rho=1e-3, verbose=True, normalize=True, norm="l2", init_coupling=True):
        if normalize:
            self.normalize(norm=norm)
        self.init_marginals()
        if init_coupling:
            self.init_coupling()
        self.init_distances(k)
        self.find_correspondences(e=e, balanced=balanced,  rho=rho, verbose=verbose)
        X_aligned, y_aligned = self.barycentric_projection()
        return X_aligned, y_aligned

    def search_scot(self, ks, es, rhos = [1], balanced = True, all_values = False, normalize=True, norm="l2"):
        '''
        Performs a hyperparameter sweep for given values of k and epsilon
        Default: return the parameters corresponding to the lowest GW distance
        (Optional): return all k, epsilon, and GW values
        '''

        # initialize alignment
        if normalize:
            self.normalize(norm=norm)
        self.init_marginals()
        self.init_coupling()
        X_aligned, y_aligned = self.barycentric_projection()

        # store values of k, epsilon, and gw distance
        total=len(es)*len(ks)*len(rhos)
        k_plot=np.zeros(total)
        e_plot=np.zeros(total)
        rho_plot=np.zeros(total)
        g_plot=np.zeros((self.num_y, total))
        gmin = [1 for i in range(self.num_y)]
        e_best =  [es[0] for i in range(self.num_y)]
        k_best =  [ks[0] for i in range(self.num_y)]
        rho_best =  [rhos[0] for i in range(self.num_y)]

        counter=0

        # search in k first to reduce graph computation
        for k in ks:
            self.init_distances(k)
            for rho in rhos:
                self.init_coupling() # reinitialize coupling for new value of epsilon
                for e in es:
                    # run scot
                    self.find_correspondences(e=e, balanced=balanced, rho=rho)
                    # save values
                    if self.flag:
                        k_plot[counter] = k
                        e_plot[counter] = e
                        rho_plot[counter] = rho

                        for i in range(self.num_y):
                            g_plot[i, counter] = self.gwdist[i]

                            # save the alignment if it is lower
                            if g_plot[i, counter] < gmin[i]:
                                X_aligned, y_aligned = self.barycentric_projection()
                                gmin[i] = g_plot[i, counter]
                                e_best[i] = e
                                k_best[i] = k
                                rho_best[i] = rho
                    
                    counter = counter + 1
           
        if all_values:
            # return alignment and all values
            return X_aligned, y_aligned, g_plot, k_plot, e_plot, rho_plot
        else:
            # return  alignment and the parameters corresponding to the lowest GW distance
            return X_aligned, y_aligned, gmin, k_best, e_best, rho_best


    def unsupervised_scot(self, balanced = True, normalize=True, norm='l2', all_values = False):
        '''
        Unsupervised hyperparameter tuning algorithm to find an alignment
        by using the GW distance as a measure of alignment
        '''

        if normalize:
            self.normalize(norm=norm)

        # use k = 20% of # sample or k = 50 if dataset is large
        n = self.X.shape[0]
        for i in range(self.num_y):
            if n > self.y[i].shape[0]:
                n = self.y[i].shape[0]
        k_start = min(n // 5, 50)

        if balanced:
            num_eps = 12
            num_k = 5
            num_rho = 1
        else:
            num_eps = 5
            num_k = 5
            num_rho = 5

        # define search space
        es = np.logspace(-1, -3, num_eps)
        if ( n > 250):
            ks = np.linspace(20, 100, num_k)
        else:
            ks = np.linspace(n//20, n//6, num_k)
        ks = ks.astype(int)
        rhos = np.logspace(-1,-3, num_rho)

        # search parameter space
        X_aligned, y_aligned, g_best, k_best, e_best, rho_best = self.search_scot(ks, 
            es, rhos, balanced=balanced, all_values=all_values, normalize=False)

        if all_values:
            return X_aligned, y_aligned, g_best, k_best, e_best, rho_best
        else:
            return X_aligned, y_aligned

