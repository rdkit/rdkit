#  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met: 
#
#     * Redistributions of source code must retain the above copyright 
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following 
#       disclaimer in the documentation and/or other materials provided 
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
#       nor the names of its contributors may be used to endorse or promote 
#       products derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#


import unittest
import os,sys

import pickle

from rdkit import rdBase
from rdkit import Chem
from rdkit.Chem import rdChemReactions, AllChem
from rdkit import Geometry
from rdkit import RDConfig
import itertools, time

test_data = [("good", '''$RXN

      ISIS     052820091627

  2  1
$MOL

  -ISIS-  05280916272D

  2  1  0  0  0  0  0  0  0  0999 V2000
   -3.2730   -7.0542    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -3.9875   -7.4667    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0
  1  2  1  0  0  0  0
V    1 halogen.bromine.aromatic
M  RGP  1   2   1
M  END
$MOL

  -ISIS-  05280916272D

  4  3  0  0  0  0  0  0  0  0999 V2000
    3.4375   -7.7917    0.0000 R#  0  0  0  0  0  0  0  0  0  2  0  0
    4.1520   -7.3792    0.0000 B   0  0  0  0  0  0  0  0  0  0  0  0
    4.1520   -6.5542    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8664   -7.7917    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0  0  0  0
  1  2  1  0  0  0  0
  2  4  1  0  0  0  0
V    2 boronicacid
M  RGP  1   1   2
M  END
$MOL

  -ISIS-  05280916272D

  2  1  0  0  0  0  0  0  0  0999 V2000
   11.2667   -7.3417    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0
   11.9811   -6.9292    0.0000 R#  0  0  0  0  0  0  0  0  0  2  0  0
  1  2  1  0  0  0  0
M  RGP  2   1   1   2   2
M  END'''),
             
("bad", '''$RXN

      ISIS     052820091627

  2  1
$MOL

  -ISIS-  05280916272D

  2  1  0  0  0  0  0  0  0  0999 V2000
   -3.2730   -7.0542    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -3.9875   -7.4667    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
V    1 halogen.bromine.aromatic
M  RGP  1   2   1
M  END
$MOL

  -ISIS-  05280916272D

  4  3  0  0  0  0  0  0  0  0999 V2000
    3.4375   -7.7917    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
    4.1520   -7.3792    0.0000 B   0  0  0  0  0  0  0  0  0  0  0  0
    4.1520   -6.5542    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8664   -7.7917    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0  0  0  0
  1  2  1  0  0  0  0
  2  4  1  0  0  0  0
V    2 boronicacid
M  RGP  1   1   2
M  END
$MOL

  -ISIS-  05280916272D

  2  1  0  0  0  0  0  0  0  0999 V2000
   11.2667   -7.3417    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
   11.9811   -6.9292    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  RGP  2   1   1   2   2
M  END'''),
# chemdraw style             
("bad", '''$RXN

      ISIS     052820091627

  2  1
$MOL

  -ISIS-  05280916272D

  2  1  0  0  0  0  0  0  0  0999 V2000
   -3.2730   -7.0542    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -3.9875   -7.4667    0.0000 R1  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
V    1 halogen.bromine.aromatic
M  END
$MOL

  -ISIS-  05280916272D

  4  3  0  0  0  0  0  0  0  0999 V2000
    3.4375   -7.7917    0.0000 R2  0  0  0  0  0  0  0  0  0  0  0  0
    4.1520   -7.3792    0.0000 B   0  0  0  0  0  0  0  0  0  0  0  0
    4.1520   -6.5542    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8664   -7.7917    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0  0  0  0
  1  2  1  0  0  0  0
  2  4  1  0  0  0  0
V    2 boronicacid
M  END
$MOL

  -ISIS-  05280916272D

  2  1  0  0  0  0  0  0  0  0999 V2000
   11.2667   -7.3417    0.0000 R1  0  0  0  0  0  0  0  0  0  0  0  0
   11.9811   -6.9292    0.0000 R2  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END'''),
("fail", '''$RXN

      ISIS     052820091627

  2  1
$MOL

  -ISIS-  05280916272D

  2  1  0  0  0  0  0  0  0  0999 V2000
   -3.2730   -7.0542    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -3.9875   -7.4667    0.0000 R1  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
V    1 halogen.bromine.aromatic
M  END
$MOL

  -ISIS-  05280916272D

  4  3  0  0  0  0  0  0  0  0999 V2000
    3.4375   -7.7917    0.0000 R3  0  0  0  0  0  0  0  0  0  0  0  0
    4.1520   -7.3792    0.0000 B   0  0  0  0  0  0  0  0  0  0  0  0
    4.1520   -6.5542    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.8664   -7.7917    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0  0  0  0
  1  2  1  0  0  0  0
  2  4  1  0  0  0  0
V    2 boronicacid
M  END
$MOL

  -ISIS-  05280916272D

  2  1  0  0  0  0  0  0  0  0999 V2000
   11.2667   -7.3417    0.0000 R1  0  0  0  0  0  0  0  0  0  0  0  0
   11.9811   -6.9292    0.0000 R2  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END'''),
]

unused_rlabel_in_product = """$RXN
bug.rxn
  ChemDraw06121709062D

  1  1
$MOL



  2  1  0  0  0  0  0  0  0  0999 V2000
    0.1604    0.3798    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1604   -0.3798    0.0000 R   0  0  0  0  0  0  0  0  0  1  0  0
  1  2  1  0        0
M  END
$MOL



  2  1  0  0  0  0  0  0  0  0999 V2000
   -1.2690   -1.3345    0.0000 R   0  0  0  0  0  0  0  0  0  1  0  0
    1.2690    1.3345    0.0000 R1  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0        0
M  END
"""
good_res = (0,0,2,1,(((0, 'halogen.bromine.aromatic'),), ((1, 'boronicacid'),)))
bad_res = (3,0,2,1,(((0, 'halogen.bromine.aromatic'),), ((1, 'boronicacid'),)))

class TestCase(unittest.TestCase) :
    def test_sanitize(self):
        for status, block in test_data:
            print("*"*44)
            rxna = AllChem.ReactionFromRxnBlock(block)
            rxnb = AllChem.ReactionFromRxnBlock(block)
            rxna.Initialize()
            res = rdChemReactions.PreprocessReaction(rxna)
            print(AllChem.ReactionToRxnBlock(rxna))
            if status == "good":
                self.assertEquals(res, good_res)
            elif status == "bad":
                self.assertEquals(res, bad_res)
            print (">"*44)
            rxnb.Initialize()
            try:
                rdChemReactions.SanitizeRxn(rxnb)
                res = rdChemReactions.PreprocessReaction(rxnb)
                print(AllChem.ReactionToRxnBlock(rxnb))
                self.assertEquals(res, good_res)
                assert not status == "fail"
            except:
                print ("$RXN Failed")
                if status == "fail":
                    continue
                raise
                
    def test_unused_rlabel_in_product(self):
        rxn = AllChem.ReactionFromRxnBlock(unused_rlabel_in_product)
        # test was for a seg fault
        rdChemReactions.SanitizeRxn(rxn)



if __name__ == '__main__':
  unittest.main()
