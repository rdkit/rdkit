## Automatically adapted for numpy.oldnumeric Jun 27, 2008 by -c

#
#  Copyright (C) 2002-2008  greg Landrum and Rational Discovery LLC
#

""" unit tests for the model and descriptor packager """
from rdkit import RDConfig
from rdkit.ML.Data import DataUtils
import unittest,os,sys
import io
from rdkit.six.moves import cPickle
from rdkit.ML.ModelPackage import Packager
from rdkit import Chem
import random

def feq(a,b,tol=1e-4):
  return abs(a-b)<=tol


class TestCase(unittest.TestCase):
  def setUp(self):
    self.dataDir =os.path.join(RDConfig.RDCodeDir,'ML/ModelPackage/test_data')
    self.testD = [
      # NOTE: the confidences here can be twitchy due to changes in descriptors:
      ('Fc1ccc(NC(=O)c2cccnc2Oc3cccc(c3)C(F)(F)F)c(F)c1',0,0.8 ),
      #(r'CN/1(=C\C=C(/C=C1)\C\2=C\C=N(C)(Cl)\C=C2)Cl',0,0.70),
      (r'NS(=O)(=O)c1cc(ccc1Cl)C2(O)NC(=O)c3ccccc32',1,0.70),
      ]

  def _verify(self,pkg,testD):
    for smi,pred,conf in testD:
      try:
        m = Chem.MolFromSmiles(smi)
      except:
        sys.stderr.write('SMILES: %s failed\n'%(smi))
      else:
        p,c = pkg.Classify(m)
        assert p==pred,'bad prediction (%d) for smiles %s'%(p,smi)
        assert feq(c,conf),'bad confidence (%f) for smiles %s'%(c,smi)
  def _verify2(self,pkg,testD):
    for smi,pred,conf in testD:
      try:
        m = Chem.MolFromSmiles(smi)
      except:
        sys.stderr.write('SMILES: %s failed\n'%(smi))
      else:
        p,c = pkg.Classify(m)
        assert p==pred,'bad prediction (%d) for smiles %s'%(p,smi)
        assert feq(c,conf),'bad confidence (%f) for smiles %s'%(c,smi)
        p,c = pkg.Classify(m)
        assert p==pred,'bad prediction (%d) for smiles %s'%(p,smi)
        assert feq(c,conf),'bad confidence (%f) for smiles %s'%(c,smi)
      
  def testBuild(self):
    """ tests building and screening a packager """
    with open(os.path.join(self.dataDir,'Jan9_build3_calc.dsc'),'r') as calcTF:
      buf = calcTF.read().replace('\r\n', '\n').encode('utf-8')
      calcTF.close()
    with io.BytesIO(buf) as calcF:
      calc = cPickle.load(calcF)
    with open(os.path.join(self.dataDir,'Jan9_build3_model.pkl'),'rb') as modelF:
      model = cPickle.load(modelF)
    pkg = Packager.ModelPackage(descCalc=calc,model=model)
    self._verify(pkg,self.testD)
  
  def testLoad(self):
    """ tests loading and screening a packager """
    with open(os.path.join(self.dataDir,'Jan9_build3_pkg.pkl'),'r') as pkgTF:
      buf = pkgTF.read().replace('\r\n', '\n').encode('utf-8')
      pkgTF.close()
    with io.BytesIO(buf) as pkgF:
      pkg = cPickle.load(pkgF)
    self._verify(pkg,self.testD)
  
  def testLoad2(self):
    """ tests loading and screening a packager 2 """
    with open(os.path.join(self.dataDir,'Jan9_build3_pkg.pkl'),'r') as pkgTF:
      buf = pkgTF.read().replace('\r\n', '\n').encode('utf-8')
      pkgTF.close()
    with io.BytesIO(buf) as pkgF:
      pkg = cPickle.load(pkgF)
    self._verify2(pkg,self.testD)
  
  def testPerm1(self):
    """ tests the descriptor remapping stuff in a packager """
    from rdkit.Chem import Descriptors
    with open(os.path.join(self.dataDir,'Jan9_build3_pkg.pkl'),'r') as pkgTF:
      buf = pkgTF.read().replace('\r\n', '\n').encode('utf-8')
      pkgTF.close()
    with io.BytesIO(buf) as pkgF:
      pkg = cPickle.load(pkgF)
    calc = pkg.GetCalculator()
    names = calc.GetDescriptorNames()
    ref = {}
    DataUtils.InitRandomNumbers((23,42))
    for smi,pred,conf in self.testD:
      for desc in names:
        fn = getattr(Descriptors,desc,lambda x:777)
        m = Chem.MolFromSmiles(smi)
        ref[desc] = fn(m)

      for i in range(5):
        perm = list(names)
        random.shuffle(perm,random=random.random)

        m = Chem.MolFromSmiles(smi)
        for desc in perm:
          fn = getattr(Descriptors,desc,lambda x:777)
          val = fn(m)
          assert feq(val,ref[desc],1e-4),'%s: %s(%s): %f!=%f'%(str(perm),
                                                               smi,
                                                               desc,
                                                               val,
                                                               ref[desc])

  def testPerm2(self):
    """ tests the descriptor remapping stuff in a packager """
    with open(os.path.join(self.dataDir,'Jan9_build3_pkg.pkl'),'r') as pkgTF:
      buf = pkgTF.read().replace('\r\n', '\n').encode('utf-8')
      pkgTF.close()
    with io.BytesIO(buf) as pkgF:
      pkg = cPickle.load(pkgF)
    calc = pkg.GetCalculator()
    names = calc.GetDescriptorNames()
    DataUtils.InitRandomNumbers((23,42))
    perm = list(names)
    random.shuffle(perm,random=random.random)
    calc.simpleList = perm
    calc.descriptorNames = perm
    pkg.Init()
    self._verify(pkg,self.testD)


  
if __name__ == '__main__':
  unittest.main()
  
  
