# Author: Yutong Zhao
# Taken from the Time Machine

import numpy as np
import jax.numpy as jnp

def embed(D, dims):
    """
    Embed a distance matrix into Euclidean coordinates.

    Parameters
    ----------
    D: np.ndarray of shape (N,N)
        A symmetric pairwise distance matrix

    dims: int
        Number of dimensions we keep

    Returns
    -------
    Cartesian Coordinates given the distance matrix (N,dims)

    """
    # Generate Gramian Matrix (B)
    n = D.shape[0]                                                                       
    H = jnp.eye(n) - jnp.ones((n, n))/n
    B = -H.dot(D**2).dot(H)/2
    
    # perform SVD (vmap batch is broken)
    # https://github.com/google/jax/issues/1397
    # u, s, v = jnp.linalg.svd(B)

    # so we do it manually instead
    evals, evecs = jnp.linalg.eigh(B)

    s = jnp.abs(evals)
    # we don't need to multiply evecs by -1 since D is hermitian

    perm = jnp.argsort(s)[::-1]
    s = s[perm]
    u = evecs[:, perm]

    x = u.dot(jnp.diag(jnp.sqrt(s)))
    return x[:, :dims]/10

def generate_nxn(bounds_mat):
    """
    Generate a N x N pairwise distance matrix by random sampling the upper and lower constraints

    Parameters
    ----------
    bounds_mat: [N, N]
        A bounds matrix

    Returns
    -------
    np.ndarray of shape (N,N)
        Pairwise distance matrix

    """
    num_atoms = bounds_mat.shape[0]
    
    dij = np.zeros((num_atoms, num_atoms))
    # We can vectorize this
    # Sample from upper/lower bound matrix
    for i in range(num_atoms):
        for j in range(i, num_atoms):
            upper_bound = bounds_mat[i][j]
            lower_bound = bounds_mat[j][i]
            random_dij = np.random.uniform(lower_bound, upper_bound)
            dij[i][j] = random_dij
            dij[j][i] = random_dij
            # Warn if lower bound is greater than upper bound
            if lower_bound > upper_bound:
                print('WARNING: lower bound {} greater than upper bound {}'.format(lower_bound,upper_bound))
    
    return dij