# -*- coding: utf-8 -*-
"""
molvs.normalize
~~~~~~~~~~~~~~~

This module contains tools for normalizing molecules using reaction SMARTS patterns.

:copyright: Copyright 2016 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""


import logging

from rdkit import Chem
from rdkit.Chem import AllChem

from .utils import memoized_property


log = logging.getLogger(__name__)


class Normalization(object):
    """A normalization transform defined by reaction SMARTS."""

    def __init__(self, name, transform):
        """
        :param string name: A name for this Normalization
        :param string transform: Reaction SMARTS to define the transformation.
        """
        log.debug('Initializing Normalization: %s', name)
        self.name = name
        self.transform_str = transform

    @memoized_property
    def transform(self):
        log.debug('Loading Normalization transform: %s', self.name)
        return AllChem.ReactionFromSmarts(str(self.transform_str))

    def __repr__(self):
        return 'Normalization({!r}, {!r})'.format(self.name, self.transform_str)

    def __str__(self):
        return self.name


#: The default list of Normalization transforms.
NORMALIZATIONS = (
    # Opposite of #2.1 in InChI technical manual? Covered by RDKit Sanitization.
    Normalization('Nitro to N+(O-)=O',
                  '[N,P,As,Sb;X3:1](=[O,S,Se,Te:2])=[O,S,Se,Te:3]>>[*+1:1]([*-1:2])=[*:3]'),
    Normalization('Sulfone to S(=O)(=O)', '[S+2:1]([O-:2])([O-:3])>>[S+0:1](=[O-0:2])(=[O-0:3])'),
    Normalization('Pyridine oxide to n+O-', '[n:1]=[O:2]>>[n+:1][O-:2]'),
    Normalization('Azide to N=N+=N-', '[*,H:1][N:2]=[N:3]#[N:4]>>[*,H:1][N:2]=[N+:3]=[N-:4]'),
    Normalization('Diazo/azo to =N+=N-', '[*:1]=[N:2]#[N:3]>>[*:1]=[N+:2]=[N-:3]'),
    Normalization('Sulfoxide to -S+(O-)-',
                  '[!O:1][S+0;X3:2](=[O:3])[!O:4]>>[*:1][S+1:2]([O-:3])[*:4]'),
    # Equivalent to #1.5 in InChI technical manual
    Normalization('Phosphate to P(O-)=O',
                  '[O,S,Se,Te;-1:1][P+;D4:2][O,S,Se,Te;-1:3]>>[*+0:1]=[P+0;D5:2][*-1:3]'),
    # Equivalent to #1.8 in InChI technical manual
    Normalization('C/S+N to C/S=N+', '[C,S;X3+1:1]([NX3:2])[NX3!H0:3]>>[*+0:1]([N:2])=[N+:3]'),
    # Equivalent to #1.8 in InChI technical manual
    Normalization('P+N to P=N+', '[P;X4+1:1]([NX3:2])[NX3!H0:3]>>[*+0:1]([N:2])=[N+:3]'),
    Normalization('Normalize hydrazine-diazonium',
                  '[CX4:1][NX3H:2]-[NX3H:3][CX4:4][NX2+:5]#[NX1:6]>>[CX4:1][NH0:2]=[NH+:3][C:4][N+0:5]=[NH:6]'),
    # Equivalent to #1.3 in InChI technical manual
    Normalization('Recombine 1,3-separated charges',
                  '[N,P,As,Sb,O,S,Se,Te;-1:1]-[A+0:2]=[N,P,As,Sb,O,S,Se,Te;+1:3]>>[*-0:1]=[*:2]-[*+0:3]'),
    Normalization('Recombine 1,3-separated charges',
                  '[n,o,p,s;-1:1]:[a:2]=[N,O,P,S;+1:3]>>[*-0:1]:[*:2]-[*+0:3]'),
    Normalization('Recombine 1,3-separated charges',
                  '[N,O,P,S;-1:1]-[a:2]:[n,o,p,s;+1:3]>>[*-0:1]=[*:2]:[*+0:3]'),
    Normalization('Recombine 1,5-separated charges',
                  '[N,P,As,Sb,O,S,Se,Te;-1:1]-[A+0:2]=[A:3]-[A:4]=[N,P,As,Sb,O,S,Se,Te;+1:5]>>[*-0:1]=[*:2]-[*:3]=[*:4]-[*+0:5]'),
    Normalization('Recombine 1,5-separated charges',
                  '[n,o,p,s;-1:1]:[a:2]:[a:3]:[c:4]=[N,O,P,S;+1:5]>>[*-0:1]:[*:2]:[*:3]:[c:4]-[*+0:5]'),
    Normalization('Recombine 1,5-separated charges',
                  '[N,O,P,S;-1:1]-[c:2]:[a:3]:[a:4]:[n,o,p,s;+1:5]>>[*-0:1]=[c:2]:[*:3]:[*:4]:[*+0:5]'),
    # Conjugated cation rules taken from Francis Atkinson's standardiser. Those that can reduce aromaticity aren't included
    Normalization('Normalize 1,3 conjugated cation',
                  '[N,O;+0!H0:1]-[A:2]=[N!$(*[O-]),O;+1H0:3]>>[*+1:1]=[*:2]-[*+0:3]'),
    Normalization('Normalize 1,3 conjugated cation',
                  '[n;+0!H0:1]:[c:2]=[N!$(*[O-]),O;+1H0:3]>>[*+1:1]:[*:2]-[*+0:3]'),
    #Normalization('Normalize 1,3 conjugated cation', '[N,O;+0!H0:1]-[c:2]:[n!$(*[O-]),o;+1H0:3]>>[*+1:1]=[*:2]:[*+0:3]'),
    Normalization('Normalize 1,5 conjugated cation',
                  '[N,O;+0!H0:1]-[A:2]=[A:3]-[A:4]=[N!$(*[O-]),O;+1H0:5]>>[*+1:1]=[*:2]-[*:3]=[*:4]-[*+0:5]'),
    Normalization('Normalize 1,5 conjugated cation',
                  '[n;+0!H0:1]:[a:2]:[a:3]:[c:4]=[N!$(*[O-]),O;+1H0:5]>>[n+1:1]:[*:2]:[*:3]:[*:4]-[*+0:5]'),
    # Normalization('Normalize 1,5 conjugated cation', '[N,O;+0!H0:1]-[c:2]:[a:3]:[a:4]:[n!$(*[O-]),o;+1H0:5]>>[*+1:1]=[c:2]:[*:3]:[*:4]:[*+0:5]'),
    # Normalization('Normalize 1,5 conjugated cation', '[n;+0!H0:1]1:[a:2]:[a:3]:[a:4]:[n!$(*[O-]);+1H0:5]1>>[n+1:1]1:[*:2]:[*:3]:[*:4]:[n+0:5]1'),
    # Normalization('Normalize 1,5 conjugated cation', '[n;+0!H0:1]:[a:2]:[a:3]:[a:4]:[n!$(*[O-]);+1H0:5]>>[n+1:1]:[*:2]:[*:3]:[*:4]:[n+0:5]'),
    # Equivalent to #1.6 in InChI technical manual. RDKit Sanitization handles this for perchlorate.
    Normalization('Charge normalization', '[F,Cl,Br,I,At;-1:1]=[O:2]>>[*-0:1][O-:2]'),
    Normalization('Charge recombination', '[N,P,As,Sb;-1:1]=[C+;v3:2]>>[*+0:1]#[C+0:2]'),
)
# InChI technical manual has many additional rules that cover situations that are disallowed by RDKit


#: The default value for the maximum number of times to attempt to apply the series of normalizations.
MAX_RESTARTS = 200


class Normalizer(object):
    """A class for applying Normalization transforms.

    This class is typically used to apply a series of Normalization transforms to correct functional groups and
    recombine charges. Each transform is repeatedly applied until no further changes occur.
    """

    def __init__(self, normalizations=NORMALIZATIONS, max_restarts=MAX_RESTARTS):
        """Initialize a Normalizer with an optional custom list of :class:`~molvs.normalize.Normalization` transforms.

        :param normalizations: A list of  :class:`~molvs.normalize.Normalization` transforms to apply.
        :param int max_restarts: The maximum number of times to attempt to apply the series of normalizations (default
                                 200).
        """
        log.debug('Initializing Normalizer')
        self.normalizations = normalizations
        self.max_restarts = max_restarts

    def __call__(self, mol):
        """Calling a Normalizer instance like a function is the same as calling its normalize(mol) method."""
        return self.normalize(mol)

    def normalize(self, mol):
        """Apply a series of Normalization transforms to correct functional groups and recombine charges.

        A series of transforms are applied to the molecule. For each Normalization, the transform is applied repeatedly
        until no further changes occur. If any changes occurred, we go back and start from the first Normalization
        again, in case the changes mean an earlier transform is now applicable. The molecule is returned once the entire
        series of Normalizations cause no further changes or if max_restarts (default 200) is reached.

        :param mol: The molecule to normalize.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :return: The normalized fragment.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        log.debug('Running Normalizer')
        # Normalize each fragment separately to get around quirky RunReactants behaviour
        fragments = []
        for fragment in Chem.GetMolFrags(mol, asMols=True):
            fragments.append(self._normalize_fragment(fragment))
        # Join normalized fragments into a single molecule again
        outmol = fragments.pop()
        for fragment in fragments:
            outmol = Chem.CombineMols(outmol, fragment)
        Chem.SanitizeMol(outmol)
        return outmol

    def _normalize_fragment(self, mol):
        for n in range(self.max_restarts):
            # Iterate through Normalization transforms and apply each in order
            for normalization in self.normalizations:
                product = self._apply_transform(mol, normalization.transform)
                if product:
                    # If transform changed mol, go back to first rule and apply each again
                    log.info('Rule applied: %s', normalization.name)
                    mol = product
                    break
            else:
                # For loop finishes normally, all applicable transforms have been applied
                return mol
        # If we're still going after max_restarts (default 200), stop and warn, but still return the mol
        log.warning('Gave up normalization after %s restarts', self.max_restarts)
        return mol

    def _apply_transform(self, mol, rule):
        """Repeatedly apply normalization transform to molecule until no changes occur.

        It is possible for multiple products to be produced when a rule is applied. The rule is applied repeatedly to
        each of the products, until no further changes occur or after 20 attempts. If there are multiple unique products
        after the final application, the first product (sorted alphabetically by SMILES) is chosen.
        """
        mols = [mol]
        for n in range(20):
            products = {}
            for mol in mols:
                for product in [x[0] for x in rule.RunReactants((mol,))]:
                    if Chem.SanitizeMol(product, catchErrors=True) == 0:
                        products[Chem.MolToSmiles(product, isomericSmiles=True)] = product
            if products:
                mols = [products[s] for s in sorted(products)]
            else:
                # If n == 0, the rule was not applicable and we return None
                return mols[0] if n > 0 else None
