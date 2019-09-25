# Author: Yutong Zhao

from jax.config import config; config.update("jax_enable_x64", True) # Disable if you don't care about precision

from openforcefield.typing.engines.smirnoff import ForceField
from rdkit import Chem
from rdkit.Chem import AllChem

import forcefield 
import numpy as np
import jax
import functools
import dg
import time
import optimizer

aspirin = "O=C(C)Oc1ccccc1C(=O)O"
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
n_confs = 100
AllChem.EmbedMultipleConfs(mol, numConfs=n_confs)

# 1. Serial Optimization with scipy, run one molecule at a time.

# # Quasi-newton methods are quite slow
for cidx in range(2):
    conf_angstroms = np.array(mol.GetConformer(cidx).GetPositions())
    conf_nanometers = conf_angstroms/10
    s_e = potential(conf_nanometers)
    minimized_conf = optimizer.minimize_scipy_bfgs(potential, gradient, conf_nanometers)
    f_e = potential(minimized_conf)
    N = conf_nanometers.shape[0]
    D = conf_nanometers.shape[1]
    hess = hessian(minimized_conf)[0][0].reshape((N*D, N*D))
    print("Minimized from ", s_e, "to", f_e, "kcal/mol with final eigenvalues:", np.linalg.eigh(hess)[0])


# 2. Batched optimization with jax optimizers
all_confs = []
for cidx in range(n_confs):
    conf_angstroms = np.array(mol.GetConformer(cidx).GetPositions())
    conf_nanometers = conf_angstroms/10
    all_confs.append(conf_nanometers)

all_confs = np.array(all_confs)

batched_potential = jax.jit(jax.vmap(potential))

print("Before Batched Energies", batched_potential(all_confs))

jax_minimizer = functools.partial(optimizer.minimize_jax_momentum, potential, gradient)
batched_minimizer = jax.jit(jax.vmap(jax_minimizer))
batched_geometries = batched_minimizer(all_confs)

print("After Batched Energies", batched_potential(batched_geometries))

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

dijs = []
for _ in range(n_confs):
    dijs.append(dg.generate_nxn(mol))
dijs = np.array(dijs) # BxNxN

@jax.jit
def conf_gen(dij):
    conf = dg.embed(dij) # diagonalization and embedding of the Gramian
    return jax_minimizer(conf)

batched_conf_gen = jax.vmap(conf_gen)
minimized_confs = batched_conf_gen(dijs)
print("After Batched Energies", batched_potential(minimized_confs))

start_time = time.time()
batched_conf_gen(dijs)
end_time = time.time()
# Takes ~500 milliseconds if we include the batched SVD time
print("Minimized", n_confs, "conformers in", end_time-start_time, "seconds") 
