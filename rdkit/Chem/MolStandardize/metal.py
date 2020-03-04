# -*- coding: utf-8 -*-
"""
molvs.metal
~~~~~~~~~~~

This module contains tools for disconnecting metal atoms that are defined as covalently bonded to non-metals.

:copyright: Copyright 2016 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""


import logging

from rdkit import Chem


log = logging.getLogger(__name__)


# TODO: This won't disconnect e.g. covalent [Na]Cl...


class MetalDisconnector(object):
    """Class for breaking covalent bonds between metals and organic atoms under certain conditions."""

    def __init__(self):
        log.debug('Initializing MetalDisconnector')
        # Initialize SMARTS to identify relevant substructures
        # TODO: Use atomic numbers instead of element symbols in SMARTS to allow for isotopes?
        self._metal_nof = Chem.MolFromSmarts(
            '[Li,Na,K,Rb,Cs,Fr,Be,Mg,Ca,Sr,Ba,Ra,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Al,Ga,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi]~[N,O,F]')
        self._metal_non = Chem.MolFromSmarts(
            '[Al,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,Hf,Ta,W,Re,Os,Ir,Pt,Au]~[B,C,Si,P,As,Sb,S,Se,Te,Cl,Br,I,At]')

    def __call__(self, mol):
        """Calling a MetalDisconnector instance like a function is the same as calling its disconnect(mol) method."""
        return self.disconnect(mol)

    def disconnect(self, mol):
        """Break covalent bonds between metals and organic atoms under certain conditions.

        The algorithm works as follows:

        - Disconnect N, O, F from any metal.
        - Disconnect other non-metals from transition metals + Al (but not Hg, Ga, Ge, In, Sn, As, Tl, Pb, Bi, Po).
        - For every bond broken, adjust the charges of the begin and end atoms accordingly.

        :param mol: The input molecule.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :return: The molecule with metals disconnected.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        log.debug('Running MetalDisconnector')
        # Remove bonds that match SMARTS
        for smarts in [self._metal_nof, self._metal_non]:
            pairs = mol.GetSubstructMatches(smarts)
            rwmol = Chem.RWMol(mol)
            orders = []
            for i, j in pairs:
                # TODO: Could get the valence contributions of the bond instead of GetBondTypeAsDouble?
                orders.append(int(mol.GetBondBetweenAtoms(i, j).GetBondTypeAsDouble()))
                rwmol.RemoveBond(i, j)
            # Adjust neighbouring charges accordingly
            mol = rwmol.GetMol()
            for n, (i, j) in enumerate(pairs):
                chg = orders[n]
                atom1 = mol.GetAtomWithIdx(i)
                atom1.SetFormalCharge(atom1.GetFormalCharge() + chg)
                atom2 = mol.GetAtomWithIdx(j)
                atom2.SetFormalCharge(atom2.GetFormalCharge() - chg)
                log.info('Removed covalent bond between %s and %s',
                         atom1.GetSymbol(), atom2.GetSymbol())
        Chem.SanitizeMol(mol)
        return mol
