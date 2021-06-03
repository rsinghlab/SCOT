import torch
import numpy as np


def euclid_dist(x, y):
    return (x[:, None, :] - y[None, :, :]).norm(p=2, dim=2)


def log_ent(x):
    return x * (x + 1e-10).log() - x + 1


def quad_kl_div(pi, ref):
    return torch.sum(pi) * torch.sum(ref * log_ent(pi / ref)) \
           + (torch.sum(pi) - torch.sum(ref)) ** 2


def l2_distortion(pi, Cx, Cy):
    mu, nu = torch.sum(pi, dim=1), torch.sum(pi, dim=0)
    A = torch.einsum('ij,i,j', Cx ** 2, mu, mu)
    B = torch.einsum('ij,i,j', Cy ** 2, nu, nu)
    C = torch.sum(torch.einsum('ij,jl->il', Cx, pi) * torch.einsum('ij,jl->il', pi, Cy))
    return A + B - 2 * C


def grad_l2_distortion(pi, Cx, Cy):
    mu, nu = torch.sum(pi, dim=1), torch.sum(pi, dim=0)
    A = torch.einsum('ij,j->i', Cx ** 2, mu)
    B = torch.einsum('kl,l->k', Cy ** 2, nu)
    C = torch.einsum('ij,kj->ik', Cx, torch.einsum('kl,jl->kj', Cy, pi))
    return 2 * (A[:, None] + B[None, :] - 2 * C)


def gw_cost(pi, a, Cx, b, Cy, rho, blur):
    if blur == 0.:
        return l2_distortion(pi, Cx, Cy) + rho * quad_kl_div(torch.sum(pi, dim=1), a) \
               + rho * quad_kl_div(torch.sum(pi, dim=0), b)
    else:
        return l2_distortion(pi, Cx, Cy) + rho * quad_kl_div(torch.sum(pi, dim=1), a) \
               + rho * quad_kl_div(torch.sum(pi, dim=0), b) + blur * quad_kl_div(pi, a[:, None] * b[None, :])

def wfr_distortion(Cx, Cy):
    dist = np.abs(Cx[:,:,None,None] - Cy[None,None,:,:])
    dist = np.minimum(dist, 0.5 * np.pi)
    dist = np.log(np.cos(dist))
    dist = np.float32(dist)
    dist[np.isnan(dist)] = - float("Inf")
    dist = np.nan_to_num(dist)
    return - 2 * dist

def wfr_grad_distortion(pi, Cx, Cy):
    dist = wfr_distortion(Cx, Cy)
    return np.einsum('ijkl,jk->il', dist, pi)
