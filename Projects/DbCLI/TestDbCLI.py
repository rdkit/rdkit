# $Id$
#
#  Copyright (C) 2007  greg Landrum
#
#   @@ All Rights Reserved  @@
#
import unittest,subprocess,os
from rdkit import RDConfig
from rdkit.Dbase.DbConnection import DbConnect

class TestCase(unittest.TestCase):
  def test1Create(self):
    p = subprocess.Popen(('python', 'CreateDb.py','--dbDir=testData/bzr','--molFormat=smiles',
                          'testData/bzr.smi'))
    res=p.wait()
    self.failIf(res)
    p=None

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
                          '--doGobbi2D',
                          'testData/bzr.sdf'))
    res=p.wait()
    self.failIf(res)
    p=None

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
    p=None

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
      self.failUnless(nbrs[lbl]=='1.000',nbrs[lbl])
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
    p=None

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
    p=None

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
    p=None

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
    p=None

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
    p=None

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnlessEqual(len(lines),49)
    os.unlink('testData/bzr/search.out')
    
    if os.path.exists('/dev/null'):
      p = subprocess.Popen(('python', 'SearchDb.py','--dbDir=testData/bzr',
                            '--outF=/dev/null',
                            '--smilesOut=testData/bzr/search.out',
                            '--smarts=cncncc',))
    else:
      p = subprocess.Popen(('python', 'SearchDb.py','--dbDir=testData/bzr',
                            '--outF=testData/crud.out',
                            '--smilesOut=testData/bzr/search.out',
                            '--smarts=cncncc',))
    res=p.wait()
    self.failIf(res)
    p=None

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnlessEqual(len(lines),49)
    os.unlink('testData/bzr/search.out')
    if os.path.exists('testData/crud.out'):
      os.unlink('testData/crud.out')
    
    p = subprocess.Popen(('python', 'SearchDb.py','--dbDir=testData/bzr',
                          '--outF=testData/bzr/search.out','--negate','--smarts=cncncc',))

    res=p.wait()
    self.failIf(res)
    p=None

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnlessEqual(len(lines),114)
    os.unlink('testData/bzr/search.out')
    
  def test2_6SearchBoth(self):
    p = subprocess.Popen(('python', 'SearchDb.py','--dbDir=testData/bzr',
                          '--outF=testData/bzr/search.out','--query=activity<6.5','--smarts=cncncc'))

    res=p.wait()
    self.failIf(res)
    p=None

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnlessEqual(len(lines),5)
    os.unlink('testData/bzr/search.out')
    
    p = subprocess.Popen(('python', 'SearchDb.py','--dbDir=testData/bzr',
                          '--outF=testData/bzr/search.out','--query=activity<6.5',
                          '--smarts=cncncc','--negate'))

    res=p.wait()
    self.failIf(res)
    p=None

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnlessEqual(len(lines),25)
    os.unlink('testData/bzr/search.out')

  def test2_7SearchGobbi(self):
    self.failUnless(os.path.exists('testData/bzr/Compounds.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/AtomPairs.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/Descriptors.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/Fingerprints.sqlt'))
    
    p = subprocess.Popen(('python', 'SearchDb.py','--dbDir=testData/bzr','--molFormat=sdf','--topN=5',
                          '--outF=testData/bzr/search.out','--similarityType=Gobbi2D',
                          'testData/bzr.sdf'))
    res=p.wait()
    self.failIf(res)
    p=None

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
    self.failUnlessEqual(splitLs[0][0],'Adinazolam')
    self.failUnlessEqual(splitLs[0][3],'alpha-hydroxytriazolam')
    self.failUnlessEqual(splitLs[0][4],'0.631')
    os.unlink('testData/bzr/search.out')

    
  def test4CreateOptions(self):
    if os.path.exists('testData/bzr/Compounds.sqlt'):
      os.unlink('testData/bzr/Compounds.sqlt')
    if os.path.exists('testData/bzr/AtomPairs.sqlt'):
      os.unlink('testData/bzr/AtomPairs.sqlt')
    if os.path.exists('testData/bzr/Descriptors.sqlt'):
      os.unlink('testData/bzr/Descriptors.sqlt')
    if os.path.exists('testData/bzr/Fingerprints.sqlt'):
      os.unlink('testData/bzr/Fingerprints.sqlt')
    
    p = subprocess.Popen(('python', 'CreateDb.py','--dbDir=testData/bzr','--molFormat=smiles',
                          '--noExtras','--noSmiles',

                          'testData/bzr.smi'))
    res=p.wait()
    self.failIf(res)
    p=None

    self.failUnless(os.path.exists('testData/bzr/Compounds.sqlt'))
    self.failIf(os.path.exists('testData/bzr/AtomPairs.sqlt'))
    self.failIf(os.path.exists('testData/bzr/Descriptors.sqlt'))
    self.failIf(os.path.exists('testData/bzr/Fingerprints.sqlt'))
    
    conn = DbConnect('testData/bzr/Compounds.sqlt')
    d = conn.GetData('molecules',fields='count(*)')
    self.failUnlessEqual(d[0][0],10)
    d = conn.GetData('molecules',fields='*')
    self.failUnlessEqual(len(d),10)
    cns = [x.lower() for x in d.GetColumnNames()]
    self.failIf('smiles' in cns)

    conn=None
    d=None
    
    if os.path.exists('testData/bzr/Compounds.sqlt'):
      os.unlink('testData/bzr/Compounds.sqlt')
    if os.path.exists('testData/bzr/AtomPairs.sqlt'):
      os.unlink('testData/bzr/AtomPairs.sqlt')
    if os.path.exists('testData/bzr/Descriptors.sqlt'):
      os.unlink('testData/bzr/Descriptors.sqlt')
    if os.path.exists('testData/bzr/Fingerprints.sqlt'):
      os.unlink('testData/bzr/Fingerprints.sqlt')
    
    p = subprocess.Popen(('python', 'CreateDb.py','--dbDir=testData/bzr','--molFormat=smiles',
                          '--noSmiles','--noFingerprints','--noLayeredFps','--noMorganFps','--noPairs','--noDescriptors',
                          'testData/bzr.smi'))
    res=p.wait()
    self.failIf(res)
    p=None

    self.failUnless(os.path.exists('testData/bzr/Compounds.sqlt'))
    self.failIf(os.path.exists('testData/bzr/AtomPairs.sqlt'))
    self.failIf(os.path.exists('testData/bzr/Descriptors.sqlt'))
    self.failIf(os.path.exists('testData/bzr/Fingerprints.sqlt'))
    
    conn = DbConnect('testData/bzr/Compounds.sqlt')
    d = conn.GetData('molecules',fields='count(*)')
    self.failUnless(d[0][0]==10)
    d = conn.GetData('molecules',fields='*')
    self.failUnless(len(d)==10)
    cns = [x.lower() for x in d.GetColumnNames()]
    self.failIf('smiles' in cns)

    p = subprocess.Popen(('python', 'CreateDb.py','--dbDir=testData/bzr','--molFormat=smiles',
                          '--noProps','--noFingerprints','--noLayeredFps','--noMorganFps','--noPairs','--noDescriptors',
                          'testData/bzr.smi'))
    res=p.wait()
    self.failIf(res)
    p=None

    self.failUnless(os.path.exists('testData/bzr/Compounds.sqlt'))
    self.failIf(os.path.exists('testData/bzr/AtomPairs.sqlt'))
    self.failIf(os.path.exists('testData/bzr/Descriptors.sqlt'))
    self.failIf(os.path.exists('testData/bzr/Fingerprints.sqlt'))
    
    conn = DbConnect('testData/bzr/Compounds.sqlt')
    d = conn.GetData('molecules',fields='count(*)')
    self.failUnlessEqual(d[0][0],10)
    d = conn.GetData('molecules',fields='*')
    self.failUnlessEqual(len(d),10)
    cns = [x.lower() for x in d.GetColumnNames()]
    self.failUnless('smiles' in cns)

    p = subprocess.Popen(('python', 'CreateDb.py','--dbDir=testData/bzr','--molFormat=smiles',
                          '--noFingerprints','--noLayeredFps','--noMorganFps','--noPairs','--noDescriptors',
                          '--maxRowsCached=4',
                          'testData/bzr.smi'))
    res=p.wait()
    self.failIf(res)
    p=None

    self.failUnless(os.path.exists('testData/bzr/Compounds.sqlt'))
    self.failIf(os.path.exists('testData/bzr/AtomPairs.sqlt'))
    self.failIf(os.path.exists('testData/bzr/Descriptors.sqlt'))
    self.failIf(os.path.exists('testData/bzr/Fingerprints.sqlt'))
    
    conn = DbConnect('testData/bzr/Compounds.sqlt')
    d = conn.GetData('molecules',fields='count(*)')
    self.failUnlessEqual(d[0][0],10)
    d = conn.GetData('molecules',fields='*')
    self.failUnlessEqual(len(d),10)
    cns = [x.lower() for x in d.GetColumnNames()]
    self.failUnless('smiles' in cns)

    p = subprocess.Popen(('python', 'CreateDb.py','--dbDir=testData/bzr','--molFormat=smiles',
                          '--noFingerprints',
                          '--noPairs','--noDescriptors',
                          '--maxRowsCached=4',
                          'testData/bzr.smi'))
    res=p.wait()
    self.failIf(res)
    p=None

    self.failUnless(os.path.exists('testData/bzr/Compounds.sqlt'))
    self.failIf(os.path.exists('testData/bzr/AtomPairs.sqlt'))
    self.failIf(os.path.exists('testData/bzr/Descriptors.sqlt'))
    self.failUnless(os.path.exists('testData/bzr/Fingerprints.sqlt'))
    
  def test5TestBackwardsCompat(self):
    if os.path.exists('testData/bzr/Compounds.sqlt'):
      os.unlink('testData/bzr/Compounds.sqlt')
    if os.path.exists('testData/bzr/AtomPairs.sqlt'):
      os.unlink('testData/bzr/AtomPairs.sqlt')
    if os.path.exists('testData/bzr/Descriptors.sqlt'):
      os.unlink('testData/bzr/Descriptors.sqlt')
    if os.path.exists('testData/bzr/Fingerprints.sqlt'):
      os.unlink('testData/bzr/Fingerprints.sqlt')
    
    p = subprocess.Popen(('python', 'CreateDb.py','--dbDir=testData/bzr',
                          '--noFingerprints','--noDescriptors',
                          'testData/bzr.sdf'))
    res=p.wait()
    self.failIf(res)
    p=None

    conn = DbConnect('testData/bzr/AtomPairs.sqlt')
    curs = conn.GetCursor()
    curs.execute('create table tmp as select compound_id,atompairfp,torsionfp from atompairs')
    p = subprocess.Popen(('python', 'SearchDb.py','--dbDir=testData/bzr','--molFormat=sdf',
                          '--topN=5','--outF=testData/bzr/search.out','--similarityType=AtomPairs',
                          '--pairTableName=tmp',
                          'testData/bzr.sdf'))
    res=p.wait()
    self.failIf(res)
    p=None

    self.failUnless(os.path.exists('testData/bzr/search.out'))
    inF = file('testData/bzr/search.out','r')
    lines=inF.readlines()
    inF=None
    self.failUnlessEqual(len(lines),163)
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


    

if __name__ == '__main__':
  unittest.main()

