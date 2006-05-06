# $Id$
#
#  Copyright (C) 2002-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
"""unit testing code for the subgraph code

"""
import unittest
import Chem
raise NotImplementedError

class TestCase(unittest.TestCase):
  def _runTest(self,smi,count,useBO):
    m = Chem.MolFromSmiles(smi)
    ps = Chem.FindUniqueSubgraphsOfLengthN(m,count,useHs=1,
                                           useBO=useBO)
    return len(ps)
  
  def testNums(self):
    vals = [
      ('CC1CCC1',3,1,4),
      ('CC1CCC1',3,0,2),
      ('CO1CCC1',3,1,4),
      ('CO1CCC1',3,0,3),
      ('CC1OCC1',3,1,7),
      ('CC1OCC1',3,0,4),

      ('CO1CCC1',4,0,3),
      ('CO1CCC1',5,0,1),
      ('CO1CCC1',5,1,1),

      ]

    for tpl in vals:
      smi,count,useBO,tgt = tpl
      cnt = self._runTest(smi,count,useBO)
      assert cnt==tgt,'bad path count (%d) for example %s'%(cnt,str(tpl))
            

    
if __name__ == '__main__':
  unittest.main()

