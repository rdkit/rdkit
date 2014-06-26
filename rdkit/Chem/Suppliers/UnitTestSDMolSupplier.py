# $Id$
#
#  Copyright (C) 2002-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

"""unit testing code for the SD file handling stuff

"""
import unittest,sys,os
from rdkit import RDConfig
from rdkit import Chem
import tempfile
from rdkit.six.moves import cStringIO as StringIO
from rdkit.six import next
class TestCase(unittest.TestCase):
  def setUp(self):
    #print '\n%s: '%self.shortDescription(),
    self.fName = os.path.join(RDConfig.RDDataDir,'NCI','first_200.props.sdf')

  
  def _testReader(self):
    " tests reads using a file name "
    supp = Chem.SDMolSupplier(self.fName)
    for i in range(10):
      m = next(supp)
      assert m,'read %d failed'%i
      assert m.GetNumAtoms(),'no atoms in mol %d'%i
    i = 100
    m = supp[i-1]
    assert m,'read %d failed'%i
    assert m.GetNumAtoms(),'no atoms in mol %d'%i
    
    l = len(supp)
    assert l == 200,'bad supplier length: %d'%(l)

    i = 12
    m = supp[i-1]
    assert m,'back index %d failed'%i
    assert m.GetNumAtoms(),'no atoms in mol %d'%i


    try:
      m = supp[201]
    except IndexError:
      fail = 1
    else:
      fail = 0
    assert fail,'out of bound read did not fail'

  def test_Writer(self):
    " tests writes using a file name "
    with open(self.fName,'r') as inf:
      inD = inf.read()
    supp = Chem.SDMolSupplier(self.fName)
    outName = tempfile.mktemp('.sdf')
    writer = Chem.SDWriter(outName)
    m1 = next(supp)
    writer.SetProps(m1.GetPropNames())
    for m in supp:
      writer.write(m)
    writer.flush()
    writer = None
    with open(outName,'r') as inf:
      outD = inf.read()
    try:
      os.unlink(outName)
    except:
      import time
      time.sleep(1)
      try:
        os.unlink(outName)
      except:
        pass
    assert inD.count('$$$$')==outD.count('$$$$'),'bad nMols in output'

  def _testStreamRoundtrip(self):
    inD = open(self.fName).read()
    supp = Chem.SDMolSupplier(self.fName)
    outName = tempfile.mktemp('.sdf')
    writer = Chem.SDWriter(outName)
    m1 = next(supp)
    for m in supp:
      writer.write(m)
    writer.flush()
    writer = None
    outD = open(outName,'r').read()
    try:
      os.unlink(outName)
    except:
      import time
      time.sleep(1)
      try:
        os.unlink(outName)
      except:
        pass
    assert inD.count('$$$$')==outD.count('$$$$'),'bad nMols in output'
    io = StringIO(outD)
    supp = Chem.SDMolSupplier(stream=io)
    outD2 = supp.Dump()
    assert outD2.count('$$$$')==len(supp),'bad nMols in output'
    assert outD2.count('$$$$')==outD.count('$$$$'),'bad nMols in output'
    assert outD2==outD,'bad outd'

  def _testLazyDataRoundtrip(self):
    inD = open(self.fName).read()
    supp = Chem.SDMolSupplier(self.fName)
    outName = tempfile.mktemp('.sdf')
    writer = Chem.SDWriter(outName)
    m1 = next(supp)
    for m in supp:
      writer.write(m)
    writer.flush()
    writer = None
    outD = open(outName,'r').read()
    try:
      os.unlink(outName)
    except:
      import time
      time.sleep(1)
      try:
        os.unlink(outName)
      except:
        pass
    assert inD.count('$$$$')==outD.count('$$$$'),'bad nMols in output'
    supp = SDMolSupplier.LazySDMolSupplier(inD=outD)
    outD2 = supp.Dump()
    assert outD2.count('$$$$')==len(supp),'bad nMols in output'
    assert outD2.count('$$$$')==outD.count('$$$$'),'bad nMols in output'
    assert outD2==outD,'bad outd'

  def _testLazyIter(self):
    " tests lazy reads using the iterator interface "
    supp = SDMolSupplier.LazySDMolSupplier(fileN=self.fName)

    nDone = 0
    for mol in supp:
      assert mol,'read %d failed'%i
      assert mol.GetNumAtoms(),'no atoms in mol %d'%i
      nDone += 1
    assert nDone==200,'bad number of molecules: %d'%(nDone)

    l = len(supp)
    assert l == 200,'bad supplier length: %d'%(l)

    i = 12
    m = supp[i-1]
    assert m,'back index %d failed'%i
    assert m.GetNumAtoms(),'no atoms in mol %d'%i


    try:
      m = supp[201]
    except IndexError:
      fail = 1
    else:
      fail = 0
    assert fail,'out of bound read did not fail'




if __name__ == '__main__':
  unittest.main()

