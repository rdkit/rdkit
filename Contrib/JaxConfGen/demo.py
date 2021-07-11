# Author: Yutong Zhao

# from jax.config import config; config.update("jax_enable_x64", True) # Disable if you don't care about precision

from openforcefield.typing.engines.smirnoff import ForceField
from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit.Chem import rdDistGeom
import rdkit.DistanceGeometry as DG

import forcefield 
import numpy as np
import jax
import functools
import dg
import time
import optimizer

aspirin = "O=C(C)Oc1ccccc1C(=O)O"
aspirin = "CCCCN(CCCC)C(=O)c1nn(c(C)c1Cl)-c1ccc(cc1C(=O)N1CCc2ccccc2C1)C(=O)NS(=O)(=O)c1ccc2ccc(I)cc2c1"
mol = Chem.MolFromSmiles(aspirin)
off = ForceField("smirnoff99Frosst.offxml")

potential, params = forcefield.parameterize(mol, off)
potential = functools.partial(potential, params=params)

gradient = jax.grad(potential, argnums=(0,)) # dU/dx
hessian = jax.hessian(potential, argnums=(0,))

# JIT the kernels
potential = jax.jit(potential)
gradient = jax.jit(gradient)
hessian = jax.jit(hessian)

# Generate some sample geometries. In practice we probably want to use crappy DG based geometries
n_confs = 128
s_time = time.time()
AllChem.EmbedMultipleConfs(mol, numConfs=n_confs)

print("RDKit took:", time.time()-s_time, "seconds for", n_confs, "conformers")
# 1. Serial Optimization with scipy, run one molecule at a time.

# quasi-newton methods are quite slow
for cidx in range(2):
    conf_angstroms = np.array(mol.GetConformer(cidx).GetPositions())
    conf_nanometers = conf_angstroms/10
    s_e = potential(conf_nanometers)
    minimized_conf = optimizer.minimize_scipy_bfgs(potential, gradient, conf_nanometers)
    f_e = potential(minimized_conf)
    N = conf_nanometers.shape[0]
    D = conf_nanometers.shape[1]
    hess = hessian(minimized_conf)[0][0].reshape((N*D, N*D))
    # print("Minimized from ", s_e, "to", f_e, "kcal/mol with final eigenvalues:", np.linalg.eigh(hess)[0])


# 2. Batched optimization with jax optimizers
all_confs = []
for cidx in range(n_confs):
    conf_angstroms = np.array(mol.GetConformer(cidx).GetPositions())
    conf_nanometers = conf_angstroms/10
    all_confs.append(conf_nanometers)

all_confs = np.array(all_confs)

batched_potential = jax.jit(jax.vmap(potential))


res = batched_potential(all_confs)
print("Before Batched Energies mean", np.mean(res), "std", np.std(res))

jax_minimizer = functools.partial(optimizer.minimize_jax_adam, potential, gradient)
batched_minimizer = jax.jit(jax.vmap(jax_minimizer))
batched_geometries = batched_minimizer(all_confs)

res = batched_potential(batched_geometries)
print("After Batched Energies mean", np.mean(res), "std", np.std(res))

np.random.shuffle(all_confs)
# Once we have the batch minimizer, we time the execution
start_time = time.time()
batched_minimizer(all_confs)
end_time = time.time()

# Takes roughly ~360 milliseconds to minimize 100 conformers simultaneously on the JIT CPU
# in double precision
print("Minimized", n_confs, "conformers in", end_time-start_time, "seconds") 

# 3. Direct DG embedding
mol.RemoveAllConformers()
bounds_mat = rdDistGeom.GetMoleculeBoundsMatrix(mol)
DG.DoTriangleSmoothing(bounds_mat)

print("Full set conformers:", n_confs)
print("Generating distance matrices...")
dijs = []
for _ in range(n_confs):
    dijs.append(dg.generate_nxn(bounds_mat))
dijs = np.array(dijs) # BxNxN

print("done")

def conf_gen(dij, dims):
    conf = dg.embed(dij, dims) # diagonalization and embedding of the Gramian
    return jax_minimizer(conf)

batch_minimizer = jax.jit(jax.vmap(jax_minimizer))

def generate_batched_embedding(dims):

    cg_fn = functools.partial(conf_gen, dims=dims)
    batched_conf_gen_fn = jax.jit(jax.vmap(cg_fn))
    return batched_conf_gen_fn

# 3D
embed_3d = generate_batched_embedding(3)
minimized_confs = embed_3d(dijs)
res = batched_potential(minimized_confs)
print("3D After Batched Energies Median", np.median(res), "std", np.std(res))

for _ in range(10):
    print("Generating distance matrices...")
    dijs = []
    for _ in range(n_confs):
        dijs.append(dg.generate_nxn(bounds_mat))
    dijs = np.array(dijs) # BxNxN

    start_time = time.time()
    embed_3d(dijs)
    end_time = time.time()
    # Takes ~500 milliseconds if we include the batched SVD time
    print("Generated", n_confs, "conformers in", end_time-start_time, "seconds") 

# Minimize in 5D then re-minimize in 3D
embed_4d = generate_batched_embedding(4)
minimized_confs_4d = embed_4d(dijs)
unminimized_confs_3d = minimized_confs_4d[:, :, :3]
minimized_confs_3d = batch_minimizer(unminimized_confs_3d)
res = batched_potential(minimized_confs_3d)
print("4D/3D After Batched Energies Median", np.median(res), "std", np.std(res))
