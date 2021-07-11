# Author: Yutong Zhao
# Taken from the Time Machine


from openforcefield.typing.engines.smirnoff import ForceField
from rdkit.Chem import rdDistGeom
import rdkit.DistanceGeometry as DG


import forcefield 
import numpy as np
import jax
import functools
import dg
import time
import optimizer

import numpy as np
import jax.numpy as jnp

import bfgs

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


def generate_conformer(mol, n_confs, dims=3):
    """
    Generate conformers for the input molecule.

    Parameters
    ----------
    mol: rdkit.ROMol
        Input molecule

    n_confs: int
        Number of conformers we generate, sorted by increasing energy

    """
    off = ForceField("smirnoff99Frosst.offxml")
    potential, params = forcefield.parameterize(mol, off)


    potential = functools.partial(potential, params=params)
    gradient = jax.grad(potential, argnums=(0,)) # dU/dx

    # JIT the kernels
    potential = jax.jit(potential)
    gradient = jax.jit(gradient)

    # jax_minimizer = functools.partial(optimizer.minimize_jax_adam, potential, gradient)
    jax_minimizer = functools.partial(optimizer.minimize_jax_fire, potential, gradient)

    batched_minimizer = jax.jit(jax.vmap(jax_minimizer))

    bounds_mat = rdDistGeom.GetMoleculeBoundsMatrix(mol)
    DG.DoTriangleSmoothing(bounds_mat)

    def conf_gen(dij):
        conf = dg.embed(dij, dims) # diagonalization and embedding of the Gramian
        return jax_minimizer(conf)

    

    batched_conf_gen_fn = jax.jit(jax.vmap(conf_gen))

    dijs = []
    for _ in range(n_confs):
        dijs.append(dg.generate_nxn(bounds_mat))
    dijs = np.array(dijs) # BxNxN

    batched_potential = jax.jit(jax.vmap(potential))
    minimized_confs = batched_conf_gen_fn(dijs)
    energies = batched_potential(minimized_confs)
    energies = np.array(energies) # convert to numpy array so we can fancy index
    energies[np.isnan(energies)] = np.inf
    # x[x == -inf] = 0
    perm = np.argsort(energies)
    sorted_min_confs = minimized_confs[perm]
    return sorted_min_confs

