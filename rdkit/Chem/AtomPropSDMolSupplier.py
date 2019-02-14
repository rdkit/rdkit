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
import textwrap


def ApplyMolListPropsToAtoms(mol, prefix, propType):
    """
    >>> m = Chem.MolFromSmiles("COC")
    >>> m.SetProp("atom.iprop.foo","1 6 9")
    >>> ApplyMolListPropsToAtoms(m,"atom.iprop.",int)
    >>> m.GetAtomWithIdx(1).GetIntProp("foo")
    6
    >>> m = Chem.MolFromSmiles("COC")
    >>> m.SetProp("atom.iprop.foo","1  9")
    >>> ApplyMolListPropsToAtoms(m,"atom.iprop.",int)
    >>> m.GetAtomWithIdx(1).HasProp("foo")
    0
    """
    if propType == str:
        setFn = "SetProp"
    elif propType == int:
        setFn = "SetIntProp"
    elif propType == float:
        setFn = "SetDoubleProp"
    elif propType == bool:
        # we need to convert the string "0" or "1" to an int and then use SetBoolProp
        propType = int
        setFn = "SetBoolProp"
    else:
        raise ValueError('bad propType')
    pns = (pn for pn in mol.GetPropNames() if pn.find(prefix) == 0)
    for pn in pns:
        pval = mol.GetProp(pn)
        splitV = pval.replace('\n', ' ').split(' ')
        if len(splitV) != mol.GetNumAtoms():
            warnings.warn(f"property value {pn} has a value list of the wrong length. Ignoring it")
            continue
        apn = pn.replace(prefix, '')
        for i, sv in enumerate(splitV):
            if sv != '':
                sv = propType(sv)
                setter = getattr(mol.GetAtomWithIdx(i), setFn)
                setter(apn, sv)


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
    ... >  <atom.dprop.PartialCharge>  (1) 
    ... 0.008 -0.314 0.008
    ... 
    ... >  <atom.iprop.NumHeavyNeighbors>  (1) 
    ... 2 2 2
    ... 
    ... >  <atom.prop.AtomLabel>  (1) 
    ... C1 N2 C3
    ... 
    ... >  <atom.bprop.IsCarbon>  (1) 
    ... 1 0 1
    ... 
    ... >  <atom.prop.PartiallyMissing>  (1) 
    ... one  three
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
    >>> m.GetAtomWithIdx(1).GetBoolProp('IsCarbon')
    0
    >>> m.GetAtomWithIdx(0).HasProp('PartiallyMissing')
    1
    >>> m.GetAtomWithIdx(1).HasProp('PartiallyMissing')
    0

    >>> m2 = suppl[0]
    >>> m2.GetNumAtoms()
    3
    >>> m2.GetAtomWithIdx(0).GetDoubleProp('PartialCharge')
    0.008

    """

    def _updateItem(self, res):
        if res is not None:
            ApplyMolListPropsToAtoms(res, 'atom.prop.', str)
            ApplyMolListPropsToAtoms(res, 'atom.iprop.', int)
            ApplyMolListPropsToAtoms(res, 'atom.dprop.', float)
            ApplyMolListPropsToAtoms(res, 'atom.bprop.', bool)
        return res

    def __next__(self):
        res = Chem.SDMolSupplier.__next__(self)
        return self._updateItem(res)

    def __getitem__(self, which):
        res = Chem.SDMolSupplier.__getitem__(self, which)
        return self._updateItem(res)


def CreateAtomListProp(mol, storeName, propVals=None, propName=None, propType=str,
                       defaultVal=None, lineLen=200):
    '''

    >>> m = Chem.MolFromSmiles('C1OC1C')
    >>> for i,at in enumerate(m.GetAtoms()):
    ...   at.SetProp('textvalue',f'atom_{i+1}')
    >>> CreateAtomListProp(m,'atomtextvalue',propName='textvalue')
    >>> for i,at in enumerate(m.GetAtoms()):
    ...   at.SetBoolProp('IsCarbon',at.GetAtomicNum()==6)
    >>> CreateAtomListProp(m,'IsCarbon',propName='IsCarbon', propType=bool)
    >>> m.GetPropsAsDict()
    {'atom.prop.atomtextvalue': 'atom_1 atom_2 atom_3 atom_4', 'atom.bprop.IsCarbon': '1 0 1 1'}
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
    >>> newmol.GetAtomWithIdx(0).GetBoolProp("IsCarbon")
    1


    by default missing property values produce a KeyError:
    >>> m = Chem.MolFromSmiles('C1OC1')
    >>> m.GetAtomWithIdx(0).SetProp('foo','bar')
    >>> m.GetAtomWithIdx(2).SetProp('foo','baz')
    >>> CreateAtomListProp(m,'foo',propName='foo')
    Traceback (most recent call last):
        ...
    KeyError: 'foo'

    But we can provide a default value to handle that:
    >>> CreateAtomListProp(m,'foo',propName='foo',defaultVal='n/a')
    >>> m.GetPropsAsDict()
    {'atom.prop.foo': 'bar n/a baz'}

    we can also leave missing values blank:
    >>> CreateAtomListProp(m,'foo',propName='foo',defaultVal='')
    >>> m.GetPropsAsDict()
    {'atom.prop.foo': 'bar  baz'}

    we deal properly with long lines:
    >>> m = Chem.MolFromSmiles('C1CC1'*10)
    >>> CreateAtomListProp(m,'long',propVals=[f'atom-{x}' for x in range(m.GetNumAtoms())])
    >>> txt = m.GetProp('atom.prop.long')
    >>> len(txt) > 200   # this is the max length for an SD file
    True
    >>> txt.count('\\n')
    1
    >>> ApplyMolListPropsToAtoms(m,"atom.prop.",str)
    >>> len([x.GetProp('long') for x in m.GetAtoms()])
    30

    '''
    if propVals is None and propName is None:
        raise ValueError("must provide at least propName or propVals")
    if propType == str:
        finalPropName = f"atom.prop.{storeName}"
    elif propType == int:
        finalPropName = f"atom.iprop.{storeName}"
    elif propType == float:
        finalPropName = f"atom.dprop.{storeName}"
    elif propType == bool:
        finalPropName = f"atom.bprop.{storeName}"
    else:
        raise ValueError('bad propType')
    if propName:
        # we are eventually always going to want a string value, so just use GetProp here
        if defaultVal is None:
            propVals = [x.GetProp(propName) for x in mol.GetAtoms()]
        else:
            propVals = [x.GetProp(propName) if x.HasProp(propName)
                        else str(defaultVal) for x in mol.GetAtoms()]

    if len(propVals) != mol.GetNumAtoms():
        raise ValueError("propVals should be the same length as the number of atoms")
    txt = ' '.join(propVals)
    if lineLen is not None:
        txt = textwrap.fill(txt, width=lineLen)
    mol.SetProp(finalPropName, txt)


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
