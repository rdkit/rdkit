# $Id$
#
#  Copyright (C) 2003-2006  Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" unit testing code for molecule suppliers

"""
from rdkit import RDConfig
import unittest,os
from rdkit.six import next
from rdkit.six.moves import cPickle
from rdkit import Chem


class TestCase(unittest.TestCase):
  def setUp(self):
    self._files=[]
  def tearDown(self):
    for fileN in self._files:
      try:
        os.unlink(fileN)
      except OSError:
        pass

  def test1SDSupplier(self):
    fileN = os.path.join(RDConfig.RDCodeDir,'VLib','NodeLib','test_data','NCI_aids.10.sdf')
    
    suppl = Chem.SDMolSupplier(fileN)
    ms = [x for x in suppl]
    assert len(ms)==10

    # test repeating:
    ms = [x for x in suppl]
    assert len(ms)==10
    # test reset:
    suppl.reset()
    m = next(suppl)
    assert m.GetProp('_Name')=='48'
    assert m.GetProp('NSC')=='48'
    assert m.GetProp('CAS_RN')=='15716-70-8'
    m = next(suppl)
    assert m.GetProp('_Name')=='78'
    assert m.GetProp('NSC')=='78'
    assert m.GetProp('CAS_RN')=='6290-84-2'

    suppl.reset()
    for i in range(10):
      m = next(suppl)
    try:
      m = next(suppl)
    except StopIteration:
      ok=1
    else:
      ok=0
    assert ok  
    
  def test2SmilesSupplier(self):
    fileN = os.path.join(RDConfig.RDCodeDir,'VLib','NodeLib','test_data','pgp_20.txt')
    
    suppl = Chem.SmilesMolSupplier(fileN,delimiter='\t',smilesColumn=2,
                                   nameColumn=1,titleLine=1)
    ms = [x for x in suppl]
    assert len(ms)==20

    # test repeating:
    ms = [x for x in suppl]
    assert len(ms)==20
    # test reset:
    suppl.reset()
    m = next(suppl)
    assert m.GetProp('_Name')=='ALDOSTERONE'
    assert m.GetProp('ID')=='RD-PGP-0001'
    m = next(suppl)
    assert m.GetProp('_Name')=='AMIODARONE'
    assert m.GetProp('ID')=='RD-PGP-0002'
    suppl.reset()
    for i in range(20):
      m = next(suppl)
    try:
      m = next(suppl)
    except StopIteration:
      ok=1
    else:
      ok=0
    assert ok  
    
  def test3SmilesSupplier(self):
    txt="""C1CC1,1
CC(=O)O,3
fail,4
CCOC,5
"""
    import tempfile
    fileN = tempfile.mktemp('.csv')
    with open(fileN,'w+') as f:
      f.write(txt)
    self._files.append(fileN)
    suppl = Chem.SmilesMolSupplier(fileN,delimiter=',',smilesColumn=0,
                                   nameColumn=1,titleLine=0)
    ms = [x for x in suppl]
    while ms.count(None):
      ms.remove(None)
    assert len(ms)==3

    
        
if __name__ == '__main__':
  unittest.main()

