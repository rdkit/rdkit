# Taken from the timemachine trunk
# From www.github.com/proteneer/timemachine

import numpy as onp
import jax.numpy as np
from jax.scipy.special import erf, erfc

from constants import ONE_4PI_EPS0


def delta_r(ri, rj, box=None):
    diff = ri - rj # this can be either N,N,3 or B,3
    if box is not None:
        diff -= box[2]*np.floor(np.expand_dims(diff[...,2], axis=-1)/box[2][2]+0.5)
        diff -= box[1]*np.floor(np.expand_dims(diff[...,1], axis=-1)/box[1][1]+0.5)
        diff -= box[0]*np.floor(np.expand_dims(diff[...,0], axis=-1)/box[0][0]+0.5)
    return diff

def distance(ri, rj, box=None):
    dxdydz = np.power(delta_r(ri, rj, box), 2)
    # np.linalg.norm nans but this doesn't
    dij = np.sqrt(np.sum(dxdydz, axis=-1))
    return dij


def harmonic_bond(conf, params, box, bond_idxs, param_idxs):
    """
    Compute the harmonic bond energy given a collection of molecules.

    This implements a harmonic angle potential: V(t) = k*(t - t0)^2 or V(t) = k*(cos(t)-cos(t0))^2

    Parameters:
    -----------
    conf: shape [num_atoms, 3] np.array
        atomic coordinates

    params: shape [num_params,] np.array
        unique parameters

    box: shape [3, 3] np.array
        periodic boundary vectors, if not None

    bond_idxs: [num_bonds, 2] np.array
        each element (src, dst) is a unique bond in the conformation

    param_idxs: [num_bonds, 2] np.array
        each element (k_idx, r_idx) maps into params for bond constants and ideal lengths

    """
    ci = conf[bond_idxs[:, 0]]
    cj = conf[bond_idxs[:, 1]]
    dij = distance(ci, cj, box)
    kbs = params[param_idxs[:, 0]]
    r0s = params[param_idxs[:, 1]]
    energy = np.sum(kbs/2 * np.power(dij - r0s, 2.0))
    return energy


def harmonic_angle(conf, params, box, angle_idxs, param_idxs, cos_angles=True):
    """
    Compute the harmonic bond energy given a collection of molecules.

    This implements a harmonic angle potential: V(t) = k*(t - t0)^2 or V(t) = k*(cos(t)-cos(t0))^2

    Parameters:
    -----------
    conf: shape [num_atoms, 3] np.array
        atomic coordinates

    params: shape [num_params,] np.array
        unique parameters

    box: shape [3, 3] np.array
        periodic boundary vectors, if not None

    angle_idxs: shape [num_angles, 3] np.array
        each element (a, b, c) is a unique angle in the conformation. atom b is defined
        to be the middle atom.

    param_idxs: shape [num_angles, 2] np.array
        each element (k_idx, t_idx) maps into params for angle constants and ideal angles

    cos_angles: True (default)
        if True, then this instead implements V(t) = k*(cos(t)-cos(t0))^2. This is far more
        numerically stable when the angle is pi.

    """
    ci = conf[angle_idxs[:, 0]]
    cj = conf[angle_idxs[:, 1]]
    ck = conf[angle_idxs[:, 2]]

    kas = params[param_idxs[:, 0]]
    a0s = params[param_idxs[:, 1]]

    vij = delta_r(ci, cj, box)
    vjk = delta_r(ck, cj, box)

    top = np.sum(np.multiply(vij, vjk), -1)
    bot = np.linalg.norm(vij, axis=-1)*np.linalg.norm(vjk, axis=-1)

    tb = top/bot

    # (ytz): we used the squared version so that we make this energy being strictly positive
    if cos_angles:
        energies = kas/2*np.power(tb - np.cos(a0s), 2)
    else:
        angle = np.arccos(tb)
        energies = kas/2*np.power(angle - a0s, 2)

    return np.sum(energies, -1)  # reduce over all angles


def signed_torsion_angle(ci, cj, ck, cl):
    """
    Batch compute the signed angle of a torsion angle.  The torsion angle
    between two planes should be periodic but not necessarily symmetric.

    Parameters
    ----------
    ci: shape [num_torsions, 3] np.array
        coordinates of the 1st atom in the 1-4 torsion angle

    cj: shape [num_torsions, 3] np.array
        coordinates of the 2nd atom in the 1-4 torsion angle

    ck: shape [num_torsions, 3] np.array
        coordinates of the 3rd atom in the 1-4 torsion angle

    cl: shape [num_torsions, 3] np.array
        coordinates of the 4th atom in the 1-4 torsion angle

    Returns
    -------
    shape [num_torsions,] np.array
        array of torsion angles.

    """

    # Taken from the wikipedia arctan2 implementation:
    # https://en.wikipedia.org/wiki/Dihedral_angle

    # We use an identical but numerically stable arctan2
    # implementation as opposed to the OpenMM energy function to
    # avoid asingularity when the angle is zero.

    rij = delta_r(cj, ci)
    rkj = delta_r(cj, ck)
    rkl = delta_r(cl, ck)

    n1 = np.cross(rij, rkj)
    n2 = np.cross(rkj, rkl)

    lhs = np.linalg.norm(n1, axis=-1)
    rhs = np.linalg.norm(n2, axis=-1)
    bot = lhs * rhs

    y = np.sum(np.multiply(np.cross(n1, n2), rkj/np.linalg.norm(rkj, axis=-1, keepdims=True)), axis=-1)
    x = np.sum(np.multiply(n1, n2), -1)

    return np.arctan2(y, x)


def periodic_torsion(conf, params, box, torsion_idxs, param_idxs):
    """
    Compute the periodic torsional energy.

    Parameters:
    -----------
    conf: shape [num_atoms, 3] np.array
        atomic coordinates

    params: shape [num_params,] np.array
        unique parameters

    box: shape [3, 3] np.array
        periodic boundary vectors, if not None

    torsion_idxs: shape [num_torsions, 4] np.array
        indices denoting the four atoms that define a torsion

    param_idxs: shape [num_torsions, 3] np.array
        indices into the params array denoting the force constant, phase, and period
    
    """

    conf = conf[:, :3] # this is defined only in 3d

    ci = conf[torsion_idxs[:, 0]]
    cj = conf[torsion_idxs[:, 1]]
    ck = conf[torsion_idxs[:, 2]]
    cl = conf[torsion_idxs[:, 3]]

    ks = params[param_idxs[:, 0]]
    phase = params[param_idxs[:, 1]]
    period = params[param_idxs[:, 2]]
    angle = signed_torsion_angle(ci, cj, ck, cl)
    nrg = ks*(1+np.cos(period * angle - phase))
    return np.sum(nrg, axis=-1)


def lennard_jones(conf, params, box, param_idxs, scale_matrix, cutoff=None):
    """
    Implements a non-periodic LJ612 potential using the Lorentz−Berthelot combining
    rules, where sig_ij = (sig_i + sig_j)/2 and eps_ij = sqrt(eps_i * eps_j).

    Parameters
    ----------
    conf: shape [num_atoms, 3] np.array
        atomic coordinates

    params: shape [num_params,] np.array
        unique parameters

    box: shape [3, 3] np.array
        periodic boundary vectors, if not None

    param_idxs: shape [num_atoms, 2] np.array
        each tuple (sig, eps) is used as part of the combining rules

    scale_matrix: shape [num_atoms, num_atoms] np.array
        scale mask denoting how we should scale interaction e[i,j].
        The elements should be between [0, 1]. If e[i,j] is 1 then the interaction
        is fully included, 0 implies it is discarded.

    cutoff: float
        Whether or not we apply cutoffs to the system. Any interactions
        greater than cutoff is fully discarded.
    
    """
    sig = params[param_idxs[:, 0]]
    eps = params[param_idxs[:, 1]]

    sig_i = np.expand_dims(sig, 0)
    sig_j = np.expand_dims(sig, 1)
    sig_ij = (sig_i + sig_j)/2
    sig_ij_raw = sig_ij

    eps_i = np.expand_dims(eps, 0)
    eps_j = np.expand_dims(eps, 1)
    eps_ij = scale_matrix * np.sqrt(eps_i * eps_j)

    eps_ij_raw = eps_ij

    ri = np.expand_dims(conf, 0)
    rj = np.expand_dims(conf, 1)

    dij = distance(ri, rj, box)

    if cutoff is not None:
        eps_ij = np.where(dij < cutoff, eps_ij, np.zeros_like(eps_ij))

    keep_mask = scale_matrix > 0

    # (ytz): this avoids a nan in the gradient in both jax and tensorflow
    sig_ij = np.where(keep_mask, sig_ij, np.zeros_like(sig_ij))
    eps_ij = np.where(keep_mask, eps_ij, np.zeros_like(eps_ij))

    sig2 = sig_ij/dij
    sig2 *= sig2
    sig6 = sig2*sig2*sig2

    energy = 4*eps_ij*(sig6-1.0)*sig6
    energy = np.where(keep_mask, energy, np.zeros_like(energy))

    # divide by two to deal with symmetry
    return np.sum(energy)/2


def pairwise_energy(conf, box, charges, cutoff):
    """
    Numerically stable implementation of the pairwise term:
    
    eij = qi*qj/dij

    """
    qi = np.expand_dims(charges, 0) # (1, N)
    qj = np.expand_dims(charges, 1) # (N, 1)
    qij = np.multiply(qi, qj)
    ri = np.expand_dims(conf, 0)
    rj = np.expand_dims(conf, 1)
    dij = distance(ri, rj, box)

    # (ytz): trick used to avoid nans in the diagonal due to the 1/dij term.
    keep_mask = 1 - np.eye(conf.shape[0])
    qij = np.where(keep_mask, qij, np.zeros_like(qij))
    dij = np.where(keep_mask, dij, np.zeros_like(dij))
    eij = np.where(keep_mask, qij/dij, np.zeros_like(dij)) # zero out diagonals

    if cutoff is not None:
        eij = np.where(dij > cutoff, np.zeros_like(eij), eij)

    return eij

def electrostatics(conf, params, box, param_idxs, scale_matrix, cutoff=None, alpha=None, kmax=None):
    """
    Compute the electrostatic potential: sum_ij qi*qj/dij

    Parameters
    ----------
    conf: shape [num_atoms, 3] np.array
        atomic coordinates

    params: shape [num_params,] np.array
        unique parameters

    box: shape [3, 3] np.array
        periodic boundary vectors, if not None then Ewald summation is used.

    param_idxs: shape [num_atoms, 2] np.array
        each tuple (sig, eps) is used as part of the combining rules

    scale_matrix: shape [num_atoms, num_atoms] np.array
        scale mask denoting how we should scale interaction e[i,j].
        The elements should be between [0, 1]. If e[i,j] is 1 then the interaction
        is fully included, 0 implies it is discarded.

    cutoff: float
        must be less than half the periodic boundary condition for each dim

    alpha: float
        alpha term controlling the erf adjustment

    kmax: int
        number of images by which we tile out reciprocal space.

    """
    charges = params[param_idxs]

    # if we use periodic boundary conditions, then the following three parameters
    # must be set in order for Ewald to make sense.
    if box is not None:

        # note that periodic boundary conditions are subject to the following
        # convention and constraints:
        # http://docs.openmm.org/latest/userguide/theory.html#periodic-boundary-conditions
        box_lengths = np.linalg.norm(box, axis=-1)
        assert cutoff is not None and cutoff >= 0.00
        assert alpha is not None
        assert kmax is not None

        # this is an implicit assumption in the Ewald calculation. If it were any larger
        # then there may be more than N^2 number of interactions.
        if np.any(box_lengths < 2*cutoff):
            raise ValueError("Box lengths cannot be smaller than twice the cutoff.")

        return ewald_energy(conf, box, charges, scale_matrix, cutoff, alpha, kmax)

    else:    
        # non periodic electrostatics is straightforward.
        # note that we do not support reaction field approximations.
        eij = scale_matrix*pairwise_energy(conf, box, charges, cutoff)

        return ONE_4PI_EPS0*np.sum(eij)/2


def self_energy(conf, charges, alpha):
    return np.sum(ONE_4PI_EPS0 * np.power(charges, 2) * alpha/np.sqrt(np.pi))


def ewald_energy(conf, box, charges, scale_matrix, cutoff, alpha, kmax):
    eij = pairwise_energy(conf, box, charges, cutoff)

    assert cutoff is not None

    # 1. Assume scale matrix is not used at all (no exceptions, no exclusions)
    # 1a. Direct Space
    eij_direct = eij * erfc(alpha*eij)
    eij_direct = ONE_4PI_EPS0*np.sum(eij_direct)/2

    # 1b. Reciprocal Space
    eij_recip = reciprocal_energy(conf, box, charges, alpha, kmax)

    # 2. Remove over estimated scale matrix contribution scaled by erf
    eij_offset = (1-scale_matrix) * eij
    eij_offset *= erf(alpha*eij_offset)
    eij_offset = ONE_4PI_EPS0*np.sum(eij_offset)/2

    return eij_direct + eij_recip - eij_offset - self_energy(conf, charges, alpha)

def reciprocal_energy(conf, box, charges, alpha, kmax):

    assert kmax > 0
    assert box is not None
    assert alpha > 0

    recipBoxSize = (2*np.pi)/np.diag(box)

    mg = []
    lowry = 0
    lowrz = 1

    numRx, numRy, numRz = kmax, kmax, kmax

    for rx in range(numRx):
        for ry in range(lowry, numRy):
            for rz in range(lowrz, numRz):
                mg.append([rx, ry, rz])
                lowrz = 1 - numRz
            lowry = 1 - numRy

    mg = np.array(onp.array(mg))

    # lattice vectors
    ki = np.expand_dims(recipBoxSize, axis=0) * mg # [nk, 3]
    ri = np.expand_dims(conf, axis=0) # [1, N, 3]
    rik = np.sum(np.multiply(ri, np.expand_dims(ki, axis=1)), axis=-1) # [nk, N]
    real = np.cos(rik)
    imag = np.sin(rik)
    eikr = real + 1j*imag # [nk, N]
    qi = charges +0j
    Sk = np.sum(qi*eikr, axis=-1)  # [nk]
    n2Sk = np.power(np.abs(Sk), 2)
    k2 = np.sum(np.multiply(ki, ki), axis=-1) # [nk]
    factorEwald = -1/(4*alpha*alpha)
    ak = np.exp(k2*factorEwald)/k2 # [nk]
    nrg = np.sum(ak * n2Sk)
    # the following volume calculation assumes the reduced PBC convention consistent
    # with that of OpenMM
    recipCoeff = (ONE_4PI_EPS0*4*np.pi)/(box[0][0]*box[1][1]*box[2][2]) 

    return recipCoeff * nrg