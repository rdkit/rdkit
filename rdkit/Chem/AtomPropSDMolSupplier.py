#
#  Copyright (C) 2019 Greg Landrum
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from rdkit import Chem
import warnings


class AtomPropSDMolSupplier(Chem.SDMolSupplier):
    """

    >>> sdf = '''
    ...      RDKit          2D
    ... 
    ...   3  3  0  0  0  0  0  0  0  0999 V2000
    ...     0.8660    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    ...    -0.4330    0.7500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    ...    -0.4330   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    ...   1  2  1  0
    ...   2  3  1  0
    ...   3  1  1  0
    ... M  END
    ... >  <AtomDProp_PartialCharge>  (1) 
    ... 0.008 -0.314 0.008
    ... 
    ... >  <AtomIProp_NumHeavyNeighbors>  (1) 
    ... 2 2 2
    ... 
    ... >  <AtomProp_AtomLabel>  (1) 
    ... C1 N2 C3
    ... 
    ... $$$$
    ... '''
    >>> suppl = AtomPropSDMolSupplier()
    >>> suppl.SetData(sdf)
    >>> m = next(suppl)
    >>> m.GetNumAtoms()
    3
    >>> m.GetAtomWithIdx(0).GetDoubleProp('PartialCharge')
    0.008
    >>> m.GetAtomWithIdx(1).GetIntProp('NumHeavyNeighbors')
    2
    >>> m.GetAtomWithIdx(2).GetProp('AtomLabel')
    'C3'

    """

    def _setVals(self, res, prefix, convert, setFn):
        pns = (pn for pn in res.GetPropNames() if pn.find(prefix) == 0)
        for pn in pns:
            pval = res.GetProp(pn)
            splitV = pval.split(' ')
            if len(splitV) != res.GetNumAtoms():
                warnings.warn(f"property value {pn} has a value list of the wrong length")
                continue
            apn = pn.replace(prefix, '')
            for i, sv in enumerate(splitV):
                sv = convert(sv)
                setter = getattr(res.GetAtomWithIdx(i), setFn)
                setter(apn, sv)

    def __next__(self):
        res = Chem.SDMolSupplier.__next__(self)
        if res is None:
            return res
        self._setVals(res, 'AtomProp_', str, 'SetProp')
        self._setVals(res, 'AtomIProp_', int, 'SetIntProp')
        self._setVals(res, 'AtomDProp_', float, 'SetDoubleProp')
        return res


def CreateAtomProp(mol, storeName, propVals=None, propName=None, propType=str):
    '''

    >>> m = Chem.MolFromSmiles('C1OC1C')
    >>> for i,at in enumerate(m.GetAtoms()):
    ...   at.SetProp('textvalue',f'atom_{i+1}')
    >>> CreateAtomProp(m,'atomtextvalue',propName='textvalue')
    >>> m.GetPropsAsDict()
    {'AtomProp_atomtextvalue': 'atom_1 atom_2 atom_3 atom_4'}
    >>> from io import StringIO
    >>> sio = StringIO()
    >>> w = Chem.SDWriter(sio)
    >>> w.write(m)
    >>> w.flush()
    >>> sdf = sio.getvalue()
    >>> suppl = AtomPropSDMolSupplier()
    >>> suppl.SetData(sdf)
    >>> newmol = next(suppl)
    >>> newmol.GetAtomWithIdx(0).GetProp("atomtextvalue")
    'atom_1'
    >>> newmol.GetAtomWithIdx(2).GetProp("atomtextvalue")
    'atom_3'

    '''
    if propVals is None and propName is None:
        raise ValueError("must provide at least propName or propVals")
    if propType == str:
        finalPropName = f"AtomProp_{storeName}"
    elif propType == int:
        finalPropName = f"AtomIProp_{storeName}"
    elif propType == float:
        finalPropName = f"AtomDProp_{storeName}"
    else:
        0
        raise ValueError('bad propType')
    if propName:
        # we are eventually always going to want a string value, so just use GetProp here
        propVals = [x.GetProp(propName) for x in mol.GetAtoms()]
    if len(propVals) != mol.GetNumAtoms():
        raise ValueError("propVals should be the same length as the number of atoms")
    mol.SetProp(finalPropName, ' '.join(propVals))


# ------------------------------------
#
#  doctest boilerplate
#


def _runDoctests(verbose=None):  # pragma: nocover
    import sys
    import doctest
    failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=verbose)
    sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
    _runDoctests()
