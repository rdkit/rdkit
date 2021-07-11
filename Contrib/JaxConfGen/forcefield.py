# Taken from the timemachine trunk
# www.github.com/proteneer/timemachine

import sys
import rdkit
from rdkit import Chem

import simtk
import functools
import numpy as np

from openforcefield.utils import toolkits
from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.topology import ValenceDict

import potentials

import jax.numpy as jnp


def to_md_units(q):
    """
    Convert a simtk Quantity to the canonical unit system.
    """
    return q.value_in_unit_system(simtk.unit.md_unit_system)


def match_bonds(mol, triples):
    bond_idxs = []
    param_idxs = []

    for smirks, k_idx, length_idx in triples:
        bond_idxs.append([bond.src, bond.dst])
        param_idxs.append([k_idx, length_idx])

    return bond_idxs, param_idxs

def parameterize(mol, forcefield):
    """
    Parameterize an RDKit molecule with a given forcefield.

    Parameters
    ----------
    mol: RDKit.Mol
        An RDKit Molecule

    forcefield: openforcefield.Forcefield
        Forcefield used to parameterize this molecule.

    Returns
    -------
    tuple nrg_fn
        A unified potential function that takes in a geometry and params, and returns a scalar potential
        in units of kJ/mol and a corresponding opaque parameter block.

    """
    # do this in a separate pass later
    global_params = []
    global_param_groups = []

    num_atoms = mol.GetNumAtoms()

    def add_param(p, p_group):
        length = len(global_params)
        global_params.append(p)
        global_param_groups.append(p_group)
        return length

    nrg_fns = []

    for handler in forcefield._parameter_handlers.items():
        handler_name, handler_params = handler
        if handler_name == 'Bonds':

            vd = ValenceDict()
            for p in handler_params.parameters:
                k_idx, l_idx = add_param(to_md_units(p.k), 0), add_param(to_md_units(p.length), 1)
                matches = toolkits.RDKitToolkitWrapper._find_smarts_matches(mol, p.smirks)

                for m in matches:
                    vd[m] = (k_idx, l_idx)

            bond_idxs = []
            bond_param_idxs = []

            for k, v in vd.items():
                bond_idxs.append(k)
                bond_param_idxs.append(v)

            partial_fn = functools.partial(
                potentials.harmonic_bond,
                bond_idxs=np.array(bond_idxs),
                param_idxs=np.array(bond_param_idxs),
                box=None,
            )

            nrg_fns.append(partial_fn)

        elif handler_name == "Angles":

            vd = ValenceDict()
            for p in handler_params.parameters:
                k_idx, a_idx = add_param(to_md_units(p.k), 2), add_param(to_md_units(p.angle), 3)
                matches = toolkits.RDKitToolkitWrapper._find_smarts_matches(mol, p.smirks)
                for m in matches:
                    vd[m] = (k_idx, a_idx)

            angle_idxs = []
            angle_param_idxs = []

            for k, v in vd.items():
                angle_idxs.append(k)
                angle_param_idxs.append(v)

            partial_fn = functools.partial(
                potentials.harmonic_angle,
                angle_idxs=np.array(angle_idxs),
                param_idxs=np.array(angle_param_idxs),
                box=None
            )

            nrg_fns.append(partial_fn)


        elif handler_name == "ImproperTorsions":
            vd = ValenceDict()

            for all_params in handler_params.parameters:

                matches = toolkits.RDKitToolkitWrapper._find_smarts_matches(mol, all_params.smirks)

                all_k_idxs = []
                all_phase_idxs = []
                all_period_idxs = []

                for k, phase, period in zip(all_params.k, all_params.phase, all_params.periodicity):

                    # (ytz): hack
                    impdivf = 3

                    k_idx, phase_idx, period_idx = add_param(to_md_units(k/impdivf), 4), add_param(to_md_units(phase), 5), add_param(period, 6),
                    all_k_idxs.append(k_idx)
                    all_phase_idxs.append(phase_idx)
                    all_period_idxs.append(period_idx)

                for m in matches:
                    t_p = []
                    for k_idx, phase_idx, period_idx in zip(all_k_idxs, all_phase_idxs, all_period_idxs):
                        t_p.append((k_idx, phase_idx, period_idx))

                    # 3-way trefoil permutation
                    others = [m[0], m[2], m[3]]
                    for p in [(others[i], others[j], others[k]) for (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]]:
                        vd[(m[1], p[0], p[1], p[2])] = t_p

            torsion_idxs = []
            torsion_param_idxs = []

            for k, vv in vd.items():
                for v in vv:
                    torsion_idxs.append(k)
                    torsion_param_idxs.append(v)

            partial_fn = functools.partial(
                potentials.periodic_torsion,
                torsion_idxs=np.array(torsion_idxs),
                param_idxs=np.array(torsion_param_idxs),
                box=None
            )

            nrg_fns.append(partial_fn)

        elif handler_name == "ProperTorsions":

            vd = ValenceDict()
            for all_params in handler_params.parameters:
                
                matches = toolkits.RDKitToolkitWrapper._find_smarts_matches(mol, all_params.smirks)

                all_k_idxs = []
                all_phase_idxs = []
                all_period_idxs = []

                for k, phase, period, idivf in zip(all_params.k, all_params.phase, all_params.periodicity, all_params.idivf):
                    k_idx, phase_idx, period_idx = add_param(to_md_units(k/idivf), 4), add_param(to_md_units(phase), 5), add_param(period, 6),
                    all_k_idxs.append(k_idx)
                    all_phase_idxs.append(phase_idx)
                    all_period_idxs.append(period_idx)

                for m in matches:
                    t_p = []
                    for k_idx, phase_idx, period_idx in zip(all_k_idxs, all_phase_idxs, all_period_idxs):
                        t_p.append((k_idx, phase_idx, period_idx))
                    vd[m] = t_p

            torsion_idxs = []
            torsion_param_idxs = []

            for k, vv in vd.items():
                for v in vv:
                    torsion_idxs.append(k)
                    torsion_param_idxs.append(v)

            partial_fn = functools.partial(
                potentials.periodic_torsion,
                torsion_idxs=np.array(torsion_idxs),
                param_idxs=np.array(torsion_param_idxs),
                box=None
            )

            nrg_fns.append(partial_fn)

        elif handler_name == "vdW":
            # lennard jones
            vd = ValenceDict()
            for param in handler_params.parameters:
                s_idx, e_idx = add_param(to_md_units(param.sigma), 8), add_param(to_md_units(param.epsilon), 9)
                matches = toolkits.RDKitToolkitWrapper._find_smarts_matches(mol, param.smirks)
                for m in matches:
                    vd[m] = (s_idx, e_idx)

                scale_matrix = np.ones(shape=(num_atoms, num_atoms), dtype=np.float64) - np.eye(num_atoms)

                # fully exclude 1-2, 1-3, 1-4 interactions
                for (src, dst) in bond_idxs:
                    scale_matrix[src][dst] = 0
                    scale_matrix[dst][src] = 0

                for (src, _, dst) in angle_idxs:
                    scale_matrix[src][dst] = 0
                    scale_matrix[dst][src] = 0

                for (src, _, _, dst) in torsion_idxs:
                    scale_matrix[src][dst] = 0
                    scale_matrix[dst][src] = 0

            lj_param_idxs = []

            for k, v in vd.items():
                lj_param_idxs.append(v)

            partial_fn = functools.partial(
                potentials.lennard_jones,
                scale_matrix=np.array(scale_matrix),
                param_idxs=np.array(lj_param_idxs),
                box=None
            )

            nrg_fns.append(partial_fn)

    def total_potential(conf, params):

        sum_e = []
        for p in nrg_fns:
            sum_e.append(p(conf, params))

        return jnp.sum(sum_e)

    return total_potential, np.array(global_params)

    # skip charges entirely and leave it to users to process
    # process charges separately
    # model = simple_charge_model()
    # vd = ValenceDict()
    # for smirks, param in model.items():

    #     param = param

    #     c_idx = add_param(param, 7)
    #     matches = toolkits.RDKitToolkitWrapper._find_smarts_matches(mol, smirks)

    #     for m in matches:
    #         vd[m] = c_idx

    #     scale_matrix = np.ones(shape=(num_atoms, num_atoms), dtype=np.float64) - np.eye(num_atoms)
    #     # fully exclude 1-2, 1-3, tbd: 1-4
    #     for (src, dst) in bond_idxs:
    #         scale_matrix[src][dst] = 0
    #         scale_matrix[dst][src] = 0

    #     for (src, _, dst) in angle_idxs:
    #         scale_matrix[src][dst] = 0
    #         scale_matrix[dst][src] = 0

    #     for (src, _, _, dst) in torsion_idxs:
    #         scale_matrix[src][dst] = 0
    #         scale_matrix[dst][src] = 0

    # charge_param_idxs = []
    # for k, v in vd.items():
    #     charge_param_idxs.append(v)

    # # print("LIGAND NET CHARGE", np.sum(np.array(global_params)[charge_param_idxs]))


    # guest_charges = np.array(global_params)[charge_param_idxs]
    # # print("LIGAND NET CHARGE", np.sum(guest_charges))
    # offsets = np.sum(guest_charges)/guest_charges.shape[0]
    # # deltas = guest_charges - new_guest_charges
    
    # for p_idx in set(charge_param_idxs):
    #     # print("ADJUSTING", p_idx)
    #     global_params[p_idx] -= offsets
    # # print("LIGAND NET CHARGE AFTER", np.sum(np.array(global_params)[charge_param_idxs]))

    # guest_charges = np.array(global_params)[charge_param_idxs]
    # # print("LIGAND NET CHARGE", guest_charges)

    # # print("SKIPPPING ")
    # nrg_fns.append((
    #    custom_ops.Electrostatics_f64,
    #    (
    #        np.array(scale_matrix, dtype=np.int32),
    #        np.array(charge_param_idxs, dtype=np.int32)
    #    )
    # ))

    # c = mol.GetConformer(0)
    # conf = np.array(c.GetPositions(), dtype=np.float64)
    # conf = conf/10 # convert to md_units

    # masses = []
    # for atom in mol.GetAtoms():
    #     masses.append(atom.GetMass())
    # masses = np.array(masses, dtype=np.float64)




# def simple_charge_model():
#     model = {
#         "[#1:1]": 0.0157,
#         "[#1:1]-[#6X4]": 0.0157,
#         "[#1:1]-[#6X4]-[#7,#8,#9,#16,#17,#35]": 0.0157,
#         "[#1:1]-[#6X4](-[#7,#8,#9,#16,#17,#35])-[#7,#8,#9,#16,#17,#35]": 0.0157,
#         "[#1:1]-[#6X4](-[#7,#8,#9,#16,#17,#35])(-[#7,#8,#9,#16,#17,#35])-[#7,#8,#9,#16,#17,#35]": 0.0157,
#         "[#1:1]-[#6X4]~[*+1,*+2]": 0.0157,
#         "[#1:1]-[#6X3]": 0.0150,
#         "[#1:1]-[#6X3]~[#7,#8,#9,#16,#17,#35]": 0.0150,
#         "[#1:1]-[#6X3](~[#7,#8,#9,#16,#17,#35])~[#7,#8,#9,#16,#17,#35]": 0.0150,
#         "[#1:1]-[#6X2]": 0.0150,
#         "[#1:1]-[#7]": -0.157,
#         "[#1:1]-[#8]": -0.2,
#         "[#1:1]-[#16]": 0.0157,
#         "[#6:1]": 0.3860,
#         "[#6X2:1]": 0.3100,
#         "[#6X4:1]": 0.3094,
#         "[#8:1]": -0.2100,
#         "[#8X2H0+0:1]": -0.1700,
#         "[#8X2H1+0:1]": -0.2104,
#         "[#7:1]": -0.200,
#         "[#16:1]": -0.2500,
#         "[#15:1]": -0.2000,
#         "[#9:1]": -0.361,
#         "[#17:1]": -0.265,
#         "[#35:1]": -0.320,
#         "[#53:1]": 0.40,
#         "[#3+1:1]": 0.0279896,
#         "[#11+1:1]": 0.0874393,
#         "[#19+1:1]": 0.1936829,
#         "[#37+1:1]": 0.3278219,
#         "[#55+1:1]": 0.4065394,
#         "[#9X0-1:1]": 0.0033640,
#         "[#17X0-1:1]": 0.0355910,
#         "[#35X0-1:1]": 0.0586554,
#         "[#53X0-1:1]": 0.0536816
#     }
#     return model

