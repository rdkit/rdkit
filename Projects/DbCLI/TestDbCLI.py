# $Id$
#
#  Copyright (C) 2007  greg Landrum
#
#   @@ All Rights Reserved  @@
#
import unittest,subprocess,os
import RDConfig
from Dbase.DbConnection import DbConnect

class TestCase(unittest.TestCase):
  def test1Create(self):
    p = subprocess.Popen(('python', 'CreateDb.py','--dbDir=testData/bzr','--molFormat=smiles',
                          'testData/bzr.smi'))
    res=p.wait()
    self.failIf(res)

    self.failUnless(os.path.exists('testData/bzr/Compounds.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/AtomPairs.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/Descriptors.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/Fingerprints.sqlt'))
    
    conn = DbConnect('testData/bzr/Compounds.sqlt')
    d = conn.GetData('molecules',fields='count(*)')
    self.failUnless(d[0][0]==10)
    
    conn = DbConnect('testData/bzr/AtomPairs.sqlt')
    d = conn.GetData('atompairs',fields='count(*)')
    self.failUnless(d[0][0]==10)
    
    conn = DbConnect('testData/bzr/Descriptors.sqlt')
    d = conn.GetData('descriptors_v1',fields='count(*)')
    self.failUnless(d[0][0]==10)
    
    conn = DbConnect('testData/bzr/Fingerprints.sqlt')
    d = conn.GetData('rdkitfps',fields='count(*)')
    self.failUnless(d[0][0]==10)

    p = subprocess.Popen(('python', 'CreateDb.py','--dbDir=testData/bzr','--molFormat=sdf',
                          'testData/bzr.sdf'))
    res=p.wait()
    self.failIf(res)

    self.failUnless(os.path.exists('testData/bzr/Compounds.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/AtomPairs.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/Descriptors.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/Fingerprints.sqlt'))
    
    conn = DbConnect('testData/bzr/Compounds.sqlt')
    d = conn.GetData('molecules',fields='count(*)')
    self.failUnless(d[0][0]==163)
    
    conn = DbConnect('testData/bzr/AtomPairs.sqlt')
    d = conn.GetData('atompairs',fields='count(*)')
    self.failUnless(d[0][0]==163)
    
    conn = DbConnect('testData/bzr/Descriptors.sqlt')
    d = conn.GetData('descriptors_v1',fields='count(*)')
    self.failUnless(d[0][0]==163)
    
    conn = DbConnect('testData/bzr/Fingerprints.sqlt')
    d = conn.GetData('rdkitfps',fields='count(*)')
    self.failUnless(d[0][0]==163)

  def test2_1SearchFPs(self):
    self.failUnless(os.path.exists('testData/bzr/Compounds.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/AtomPairs.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/Descriptors.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/Fingerprints.sqlt'))
    
    p = subprocess.Popen(('python', 'SearchDb.py','--dbDir=testData/bzr','--molFormat=sdf',
                          '--topN=5','--outF=testData/bzr/search.out','testData/bzr.sdf'))
    res=p.wait()
    self.failIf(res)

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnless(len(lines)==163)
    splitLs=[x.strip().split(',') for x in lines]
    for line in splitLs:
      lbl = line[0]
      i=1
      nbrs={}
      lastVal=1.0
      while i<len(line):
        nbrs[line[i]]=line[i+1]
        self.failUnless(float(line[i+1])<=lastVal)
        lastVal=float(line[i+1])
        i+=2
      self.failUnless(nbrs.has_key(lbl))
      self.failUnless(nbrs[lbl]=='1.000')
    os.unlink('testData/bzr/search.out')
    
  def test2_2SearchAtomPairs(self):
    self.failUnless(os.path.exists('testData/bzr/Compounds.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/AtomPairs.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/Descriptors.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/Fingerprints.sqlt'))
    
    p = subprocess.Popen(('python', 'SearchDb.py','--dbDir=testData/bzr','--molFormat=sdf',
                          '--topN=5','--outF=testData/bzr/search.out','--similarityType=AtomPairs',
                          'testData/bzr.sdf'))
    res=p.wait()
    self.failIf(res)

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnless(len(lines)==163)
    splitLs=[x.strip().split(',') for x in lines]
    for line in splitLs:
      lbl = line[0]
      i=1
      nbrs={}
      lastVal=1.0
      while i<len(line):
        nbrs[line[i]]=line[i+1]
        self.failUnless(float(line[i+1])<=lastVal)
        lastVal=float(line[i+1])
        i+=2
      self.failUnless(nbrs.has_key(lbl))
      self.failUnless(nbrs[lbl]=='1.000')
    os.unlink('testData/bzr/search.out')
    
  def test2_3SearchTorsions(self):
    self.failUnless(os.path.exists('testData/bzr/Compounds.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/AtomPairs.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/Descriptors.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/Fingerprints.sqlt'))
    
    p = subprocess.Popen(('python', 'SearchDb.py','--dbDir=testData/bzr','--molFormat=sdf','--topN=5',
                          '--outF=testData/bzr/search.out','--similarityType=TopologicalTorsions',
                          'testData/bzr.sdf'))
    res=p.wait()
    self.failIf(res)

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnless(len(lines)==163)
    splitLs=[x.strip().split(',') for x in lines]
    for line in splitLs:
      lbl = line[0]
      i=1
      nbrs={}
      lastVal=1.0
      while i<len(line):
        nbrs[line[i]]=line[i+1]
        self.failUnless(float(line[i+1])<=lastVal)
        lastVal=float(line[i+1])
        i+=2
      self.failUnless(nbrs.has_key(lbl))
      self.failUnless(nbrs[lbl]=='1.000')
    os.unlink('testData/bzr/search.out')
    

  def test2_4SearchProps(self):
    self.failUnless(os.path.exists('testData/bzr/Compounds.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/AtomPairs.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/Descriptors.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/Fingerprints.sqlt'))
    
    p = subprocess.Popen(('python', 'SearchDb.py','--dbDir=testData/bzr',
                          '--outF=testData/bzr/search.out','--query=activity<6.5'))

    res=p.wait()
    self.failIf(res)

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnless(len(lines)==30)
    os.unlink('testData/bzr/search.out')
    
    p = subprocess.Popen(('python', 'SearchDb.py','--dbDir=testData/bzr',
                          '--outF=testData/bzr/search.out','--query=activity<6.5'))

    res=p.wait()
    self.failIf(res)

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnless(len(lines)==30)
    os.unlink('testData/bzr/search.out')
    
  def test2_5SearchSmarts(self):
    p = subprocess.Popen(('python', 'SearchDb.py','--dbDir=testData/bzr',
                          '--outF=testData/bzr/search.out','--smarts=cncncc',))

    res=p.wait()
    self.failIf(res)

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnless(len(lines)==49)
    os.unlink('testData/bzr/search.out')
    
    p = subprocess.Popen(('python', 'SearchDb.py','--dbDir=testData/bzr',
                          '--outF=/dev/null',
                          '--smilesOut=testData/bzr/search.out',
                          '--smarts=cncncc',))

    res=p.wait()
    self.failIf(res)

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnless(len(lines)==49)
    os.unlink('testData/bzr/search.out')
    
    p = subprocess.Popen(('python', 'SearchDb.py','--dbDir=testData/bzr',
                          '--outF=testData/bzr/search.out','--negate','--smarts=cncncc',))

    res=p.wait()
    self.failIf(res)

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnless(len(lines)==114)
    os.unlink('testData/bzr/search.out')
    
  def test2_6SearchBoth(self):
    p = subprocess.Popen(('python', 'SearchDb.py','--dbDir=testData/bzr',
                          '--outF=testData/bzr/search.out','--query=activity<6.5','--smarts=cncncc'))

    res=p.wait()
    self.failIf(res)

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnless(len(lines)==5)
    os.unlink('testData/bzr/search.out')
    
    p = subprocess.Popen(('python', 'SearchDb.py','--dbDir=testData/bzr',
                          '--outF=testData/bzr/search.out','--query=activity<6.5',
                          '--smarts=cncncc','--negate'))

    res=p.wait()
    self.failIf(res)

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnless(len(lines)==25)
    os.unlink('testData/bzr/search.out')

  def test3_1SDSearch(self):
    self.failUnless(os.path.exists('testData/bzr/Compounds.sqlt'))

    p = subprocess.Popen(('python', 'SDSearch.py','--dbName=testData/bzr/Compounds.sqlt',
                          '--nameOut=testData/bzr/search.out'))
    res=p.wait()
    self.failIf(res)

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnless(len(lines)==163)
    os.unlink('testData/bzr/search.out')
    
  def test3_2SDSearch(self):
    self.failUnless(os.path.exists('testData/bzr/Compounds.sqlt'))

    p = subprocess.Popen(('python', 'SDSearch.py','--dbName=testData/bzr/Compounds.sqlt',
                          '--nameOut=testData/bzr/search.out', '-q activity<6.0'))
    res=p.wait()
    self.failIf(res)

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnless(len(lines)==17)
    os.unlink('testData/bzr/search.out')
    
  def test3_3SDSearch(self):
    import Chem
    self.failUnless(os.path.exists('testData/bzr/Compounds.sqlt'))

    p = subprocess.Popen(('python', 'SDSearch.py','--dbName=testData/bzr/Compounds.sqlt',
                          '--nameOut=testData/bzr/search.out', '--smarts=cncncc',
                          '--sdOut=testData/bzr/search.out.sdf'))
    res=p.wait()
    self.failIf(res)

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnless(len(lines)==49)

    suppl=Chem.SDMolSupplier('testData/bzr/search.out.sdf')
    ms = [x for x in suppl]
    self.failUnless(len(ms)==49)
    for i,m in enumerate(ms):
      self.failUnless(m.GetProp('_Name')==lines[i].strip())

    os.unlink('testData/bzr/search.out')
    os.unlink('testData/bzr/search.out.sdf')
    
  def test3_4SDSearch(self):
    import Chem
    self.failUnless(os.path.exists('testData/bzr/Compounds.sqlt'))

    p = subprocess.Popen(('python', 'SDSearch.py','--dbName=testData/bzr/Compounds.sqlt',
                          '--nameOut=testData/bzr/search.out', '--smarts=cncncc',
                          '-q activity<6.0',
                          '--sdOut=testData/bzr/search.out.sdf'))
    res=p.wait()
    self.failIf(res)

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnless(len(lines)==5)

    suppl=Chem.SDMolSupplier('testData/bzr/search.out.sdf')
    ms = [x for x in suppl]
    self.failUnless(len(ms)==5)
    for i,m in enumerate(ms):
      self.failUnless(m.GetProp('_Name')==lines[i].strip())

    os.unlink('testData/bzr/search.out')
    os.unlink('testData/bzr/search.out.sdf')
    
    
    

if __name__ == '__main__':
  unittest.main()

