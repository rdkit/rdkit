# -*- coding: utf-8 -*-
"""
molvs.standardize
~~~~~~~~~~~~~~~~~

This module contains the main :class:`~molvs.standardize.Standardizer` class that can be used to perform all
standardization tasks, as well as convenience functions like :func:`~molvs.standardize.standardize_smiles` for common
standardization tasks.

:copyright: Copyright 2016 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""


import copy
import logging

from rdkit import Chem

from .metal import MetalDisconnector
from .fragment import PREFER_ORGANIC, LargestFragmentChooser, FragmentRemover
from .normalize import NORMALIZATIONS, MAX_RESTARTS, Normalizer
from .tautomer import TAUTOMER_TRANSFORMS, TAUTOMER_SCORES, MAX_TAUTOMERS, TautomerCanonicalizer, TautomerEnumerator
from .charge import ACID_BASE_PAIRS, CHARGE_CORRECTIONS, Reionizer, Uncharger
from .utils import memoized_property


log = logging.getLogger(__name__)


class Standardizer(object):
    """The main class for performing standardization of molecules and deriving parent molecules.

    The primary usage is via the :meth:`~molvs.standardize.Standardizer.standardize` method::

        s = Standardizer()
        mol1 = Chem.MolFromSmiles('C1=CC=CC=C1')
        mol2 = s.standardize(mol1)

    There are separate methods to derive fragment, charge, tautomer, isotope and stereo parent molecules.

    """

    def __init__(self, normalizations=NORMALIZATIONS, acid_base_pairs=ACID_BASE_PAIRS,
                 charge_corrections=CHARGE_CORRECTIONS, tautomer_transforms=TAUTOMER_TRANSFORMS,
                 tautomer_scores=TAUTOMER_SCORES, max_restarts=MAX_RESTARTS, max_tautomers=MAX_TAUTOMERS,
                 prefer_organic=PREFER_ORGANIC):
        """Initialize a Standardizer with optional custom parameters.

        :param normalizations: A list of Normalizations to apply (default: :data:`~molvs.normalize.NORMALIZATIONS`).
        :param acid_base_pairs: A list of AcidBasePairs for competitive reionization (default:
                                :data:`~molvs.charge.ACID_BASE_PAIRS`).
        :param charge_corrections: A list of ChargeCorrections to apply (default:
                                :data:`~molvs.charge.CHARGE_CORRECTIONS`).
        :param tautomer_transforms: A list of TautomerTransforms to apply (default:
                                    :data:`~molvs.tautomer.TAUTOMER_TRANSFORMS`).
        :param tautomer_scores: A list of TautomerScores used to determine canonical tautomer (default:
                                :data:`~molvs.tautomer.TAUTOMER_SCORES`).
        :param max_restarts: The maximum number of times to attempt to apply the series of normalizations (default 200).
        :param max_tautomers: The maximum number of tautomers to enumerate (default 1000).
        :param prefer_organic: Whether to prioritize organic fragments when choosing fragment parent (default False).
        """
        log.debug('Initializing Standardizer')
        self.normalizations = normalizations
        self.acid_base_pairs = acid_base_pairs
        self.charge_corrections = charge_corrections
        self.tautomer_transforms = tautomer_transforms
        self.tautomer_scores = tautomer_scores
        self.max_restarts = max_restarts
        self.max_tautomers = max_tautomers
        self.prefer_organic = prefer_organic

    def __call__(self, mol):
        """Calling a Standardizer instance like a function is the same as calling its
        :meth:`~molvs.standardize.Standardizer.standardize` method."""
        return self.standardize(mol)

    def standardize(self, mol):
        """Return a standardized version the given molecule.

        The standardization process consists of the following stages: RDKit
        :rdkit:`RemoveHs <Chem.rdmolops-module.html#RemoveHs>`, RDKit
        :rdkit:`SanitizeMol <Chem.rdmolops-module.html#SanitizeMol>`, :class:`~molvs.metal.MetalDisconnector`,
        :class:`~molvs.normalize.Normalizer`, :class:`~molvs.charge.Reionizer`, RDKit
        :rdkit:`AssignStereochemistry <Chem.rdmolops-module.html#AssignStereochemistry>`.

        :param mol: The molecule to standardize.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :returns: The standardized molecule.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        mol_props = mol.GetPropsAsDict()
        mol = copy.deepcopy(mol)
        Chem.SanitizeMol(mol)
        mol = Chem.RemoveHs(mol)
        mol = self.disconnect_metals(mol)
        mol = self.normalize(mol)
        mol = self.reionize(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        for k, v in mol_props.items():
            mol.SetProp(k, str(v))
        # TODO: Check this removes symmetric stereocenters
        return mol

    def tautomer_parent(self, mol, skip_standardize=False):
        """Return the tautomer parent of a given molecule.

        :param mol: The input molecule.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :param bool skip_standardize: Set to True if mol has already been standardized.
        :returns: The tautomer parent molecule.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        if not skip_standardize:
            mol = self.standardize(mol)
        tautomer = self.canonicalize_tautomer(mol)
        tautomer = self.standardize(tautomer)
        return tautomer

    def fragment_parent(self, mol, skip_standardize=False):
        """Return the fragment parent of a given molecule.

        The fragment parent is the largest organic covalent unit in the molecule.

        :param mol: The input molecule.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :param bool skip_standardize: Set to True if mol has already been standardized.
        :returns: The fragment parent molecule.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        if not skip_standardize:
            mol = self.standardize(mol)
        # TODO: Consider applying FragmentRemover first to remove salts, solvents?
        fragment = self.largest_fragment(mol)
        return fragment

    def stereo_parent(self, mol, skip_standardize=False):
        """Return the stereo parent of a given molecule.

        The stereo parent has all stereochemistry information removed from tetrahedral centers and double bonds.

        :param mol: The input molecule.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :param bool skip_standardize: Set to True if mol has already been standardized.
        :returns: The stereo parent molecule.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        if not skip_standardize:
            mol = self.standardize(mol)
        else:
            mol = copy.deepcopy(mol)
        Chem.RemoveStereochemistry(mol)
        return mol

    def isotope_parent(self, mol, skip_standardize=False):
        """Return the isotope parent of a given molecule.

        The isotope parent has all atoms replaced with the most abundant isotope for that element.

        :param mol: The input molecule.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :param bool skip_standardize: Set to True if mol has already been standardized.
        :returns: The isotope parent molecule.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        if not skip_standardize:
            mol = self.standardize(mol)
        else:
            mol = copy.deepcopy(mol)
        # Replace isotopes with common weight
        for atom in mol.GetAtoms():
            atom.SetIsotope(0)
        return mol

    def charge_parent(self, mol, skip_standardize=False):
        """Return the charge parent of a given molecule.

        The charge parent is the uncharged version of the fragment parent.

        :param mol: The input molecule.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :param bool skip_standardize: Set to True if mol has already been standardized.
        :returns: The charge parent molecule.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        # TODO: All ionized acids and bases should be neutralised.
        if not skip_standardize:
            mol = self.standardize(mol)
        fragment = self.fragment_parent(mol, skip_standardize=True)
        if fragment:
            uncharged = self.uncharge(fragment)
            # During final standardization, the Reionizer ensures any remaining charges are in the right places
            uncharged = self.standardize(uncharged)
            return uncharged

    def super_parent(self, mol, skip_standardize=False):
        """Return the super parent of a given molecule.

        THe super parent is fragment, charge, isotope, stereochemistry and tautomer insensitive. From the input
        molecule, the largest fragment is taken. This is uncharged and then isotope and stereochemistry information is
        discarded. Finally, the canonical tautomer is determined and returned.

        :param mol: The input molecule.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :param bool skip_standardize: Set to True if mol has already been standardized.
        :returns: The super parent molecule.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        if not skip_standardize:
            mol = self.standardize(mol)
        # We don't need to get fragment parent, because the charge parent is the largest fragment
        mol = self.charge_parent(mol, skip_standardize=True)
        mol = self.isotope_parent(mol, skip_standardize=True)
        mol = self.stereo_parent(mol, skip_standardize=True)
        mol = self.tautomer_parent(mol, skip_standardize=True)
        mol = self.standardize(mol)
        return mol

    def standardize_with_parents(self, mol):
        """"""
        standardized = self.standardize(mol)
        tautomer = self.tautomer_parent(standardized, skip_standardize=True)
        super = self.super_parent(standardized, skip_standardize=True)
        # TODO: Add other parents - have optional argument to specify which are wanted
        mols = {
            'standardized': standardized,
            'tautomer_parent': tautomer,
            'super_parent': super
        }
        return mols

    # TODO: All unique tautomers
    # TODO: All unique fragments (each has to be standardized again?)

    @memoized_property
    def disconnect_metals(self):
        """
        :returns: A callable :class:`~molvs.metal.MetalDisconnector` instance.
        """
        return MetalDisconnector()

    @memoized_property
    def normalize(self):
        """
        :returns: A callable :class:`~molvs.normalize.Normalizer` instance.
        """
        return Normalizer(normalizations=self.normalizations, max_restarts=self.max_restarts)

    @memoized_property
    def reionize(self):
        """
        :returns: A callable :class:`~molvs.charge.Reionizer` instance.
        """
        return Reionizer(acid_base_pairs=self.acid_base_pairs, charge_corrections=self.charge_corrections)

    @memoized_property
    def uncharge(self):
        """
        :returns: A callable :class:`~molvs.charge.Uncharger` instance.
        """
        return Uncharger()

    @memoized_property
    def remove_fragments(self):
        """
        :returns: A callable :class:`~molvs.fragment.FragmentRemover` instance.
        """
        return FragmentRemover()

    @memoized_property
    def largest_fragment(self):
        """
        :returns: A callable :class:`~molvs.fragment.LargestFragmentChooser` instance.
        """
        return LargestFragmentChooser(prefer_organic=self.prefer_organic)

    @memoized_property
    def enumerate_tautomers(self):
        """
        :returns: A callable :class:`~molvs.tautomer.TautomerEnumerator` instance.
        """
        return TautomerEnumerator(transforms=self.tautomer_transforms, max_tautomers=self.max_tautomers)

    @memoized_property
    def canonicalize_tautomer(self):
        """
        :returns: A callable :class:`~molvs.tautomer.TautomerCanonicalizer` instance.
        """
        return TautomerCanonicalizer(transforms=self.tautomer_transforms, scores=self.tautomer_scores,
                                     max_tautomers=self.max_tautomers)


def standardize_smiles(smiles):
    """Return a standardized canonical SMILES string given a SMILES string.

    Note: This is a convenience function for quickly standardizing a single SMILES string. It is more efficient to use
    the :class:`~molvs.standardize.Standardizer` class directly when working with many molecules or when custom options
    are needed.

    :param string smiles: The SMILES for the molecule.
    :returns: The SMILES for the standardized molecule.
    :rtype: string.
    """
    # Skip sanitize as standardize does this anyway
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol = Standardizer().standardize(mol)
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def enumerate_tautomers_smiles(smiles):
    """Return a set of tautomers as SMILES strings, given a SMILES string.

    :param smiles: A SMILES string.
    :returns: A set containing SMILES strings for every possible tautomer.
    :rtype: set of strings.
    """
    # Skip sanitize as standardize does this anyway
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol = Standardizer().standardize(mol)
    tautomers = TautomerEnumerator().enumerate(mol)
    return {Chem.MolToSmiles(m, isomericSmiles=True) for m in tautomers}


def canonicalize_tautomer_smiles(smiles):
    """Return a standardized canonical tautomer SMILES string given a SMILES string.

    Note: This is a convenience function for quickly standardizing and finding the canonical tautomer for a single
    SMILES string. It is more efficient to use the :class:`~molvs.standardize.Standardizer` class directly when working
    with many molecules or when custom options are needed.

    :param string smiles: The SMILES for the molecule.
    :returns: The SMILES for the standardize canonical tautomer.
    :rtype: string.
    """
    # Skip sanitize as standardize does this anyway
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol = Standardizer().standardize(mol)
    tautomer = TautomerCanonicalizer().canonicalize(mol)
    return Chem.MolToSmiles(tautomer, isomericSmiles=True)
