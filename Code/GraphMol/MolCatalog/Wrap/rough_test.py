# $Id$
#
#  Copyright (C) 2006  Greg Landrum
#
import unittest,os,sys
from rdkit.six.moves import cPickle
from rdkit import RDConfig
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import MolCatalog

class TestCase(unittest.TestCase):
  def test1(self):
    cat = MolCatalog.CreateMolCatalog()
    es = []
    for smi in ('C1CCC1OC','C1CCC1','C'):
      m = Chem.MolFromSmiles(smi)
      entry = MolCatalog.MolCatalogEntry()
      entry.SetMol(m)
      self.assertTrue(entry.GetMol())
      eSmi = Chem.MolToSmiles(entry.GetMol())
      self.assertTrue(eSmi==Chem.MolToSmiles(m))
      entry.SetDescription(smi)
      self.assertTrue(entry.GetDescription()==smi)
      es.append(entry)

    v=cat.AddEntry(es[0])
    self.assertTrue(v==0)
    self.assertTrue(cat.GetNumEntries()==1)

    v=cat.AddEntry(es[1])
    self.assertTrue(v==1)
    self.assertTrue(cat.GetNumEntries()==2)

    v=cat.AddEntry(es[2])
    self.assertTrue(v==2)
    self.assertTrue(cat.GetNumEntries()==3)

    cat.AddEdge(0,1)
    cat.AddEdge(0,2)
    cat.AddEdge(1,2)

    d = cPickle.dumps(cat)
    es = None
    entry = None
    cat=None

    cat = cPickle.loads(d)
    self.assertTrue(cat.GetNumEntries()==3)
    cat=None

    
     

    

if __name__ == '__main__':
    unittest.main()
