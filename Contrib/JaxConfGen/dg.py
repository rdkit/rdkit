# Author: Yutong Zhao
# Taken from the Time Machine

import numpy as np
import jax.numpy as jnp
from rdkit.Chem import rdDistGeom
import rdkit.DistanceGeometry as DG

def embed(D):
    '''
    D: distance matrix (N,N)
    returns: coordinates given the distance matrix (N,3)
    '''
    # Generate Gramian Matrix (B)
    # n = len(D)
    n = D.shape[0]                                                                       
    H = jnp.eye(n) - jnp.ones((n, n))/n
    B = -H.dot(D**2).dot(H)/2
    
    # perform SVD (vmap batch is broken)
    # u, s, v = jnp.linalg.svd(B)

    # so we do it manually instead
    evals, evecs = jnp.linalg.eigh(B)

    s = jnp.abs(evals)
    # we don't need to multiply evecs by -1 since D is hermitian

    perm = jnp.argsort(s)[::-1]
    s = s[perm]
    u = evecs[:, perm]

    x = u.dot(jnp.diag(jnp.sqrt(s)))

    return x[:, :3]/10

def generate_nxn(mol):
    '''
    mol: RDKit molecule
    returns: Conformation based only on the bounded distance matrix (shape N,3)
    '''
    bounds_mat = rdDistGeom.GetMoleculeBoundsMatrix(mol)
    
    DG.DoTriangleSmoothing(bounds_mat)
    
    num_atoms = mol.GetNumAtoms()
    
    dij = np.zeros((num_atoms, num_atoms))
    
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
    # conf = cmdscale(dij)
    
    # return conf
