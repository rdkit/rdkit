# -*- coding: utf-8 -*-
"""
molvs.resonance
~~~~~~~~~~~~~~~

Resonance (mesomeric) transformations.

:copyright: Copyright 2016 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

import logging

from rdkit import Chem


log = logging.getLogger(__name__)


MAX_STRUCTURES = 1000


class ResonanceEnumerator(object):
    """Simple wrapper around RDKit ResonanceMolSupplier.

    """

    def __init__(self, kekule_all=False, allow_incomplete_octets=False, unconstrained_cations=False,
                 unconstrained_anions=False, allow_charge_separation=False, max_structures=MAX_STRUCTURES):
        """

        :param bool allow_incomplete_octets: include resonance structures whose octets are less complete than the most octet-complete structure.
        :param bool allow_charge_separation: include resonance structures featuring charge separation also when uncharged resonance structures exist.
        :param bool kekule_all: enumerate all possible degenerate Kekule resonance structures (the default is to include just one).
        :param bool unconstrained_cations: if False positively charged atoms left and right of N with an incomplete octet are acceptable only if the conjugated group has a positive total formal charge.
        :param bool unconstrained_anions: if False, negatively charged atoms left of N are acceptable only if the conjugated group has a negative total formal charge.
        :param int max_structures: Maximum number of resonance forms.
        """
        self.kekule_all = kekule_all
        self.allow_incomplete_octets = allow_incomplete_octets
        self.unconstrained_cations = unconstrained_cations
        self.unconstrained_anions = unconstrained_anions
        self.allow_charge_separation = allow_charge_separation
        self.max_structures = max_structures

    def __call__(self, mol):
        """Calling a ResonanceEnumerator instance like a function is the same as calling its enumerate(mol) method."""
        return self.enumerate(mol)

    def enumerate(self, mol):
        """Enumerate all possible resonance forms and return them as a list.

        :param mol: The input molecule.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :return: A list of all possible resonance forms of the molecule.
        :rtype: list of :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        flags = 0
        if self.kekule_all:
            flags = flags | Chem.KEKULE_ALL
        if self.allow_incomplete_octets:
            flags = flags | Chem.ALLOW_INCOMPLETE_OCTETS
        if self.allow_charge_separation:
            flags = flags | Chem.ALLOW_CHARGE_SEPARATION
        if self.unconstrained_anions:
            flags = flags | Chem.UNCONSTRAINED_ANIONS
        if self.unconstrained_cations:
            flags = flags | Chem.UNCONSTRAINED_CATIONS
        results = []
        for result in Chem.ResonanceMolSupplier(mol, flags=flags, maxStructs=self.max_structures):
            # This seems necessary? ResonanceMolSupplier only does a partial sanitization
            Chem.SanitizeMol(result)
            results.append(result)
        return results

        # Potentially interesting: getNumConjGrps(), getBondConjGrpIdx() and getAtomConjGrpIdx()


def enumerate_resonance_smiles(smiles):
    """Return a set of resonance forms as SMILES strings, given a SMILES string.

    :param smiles: A SMILES string.
    :returns: A set containing SMILES strings for every possible resonance form.
    :rtype: set of strings.
    """
    mol = Chem.MolFromSmiles(smiles)
    # Chem.SanitizeMol(mol)  # MolFromSmiles does Sanitize by default
    mesomers = ResonanceEnumerator().enumerate(mol)
    return {Chem.MolToSmiles(m, isomericSmiles=True) for m in mesomers}
