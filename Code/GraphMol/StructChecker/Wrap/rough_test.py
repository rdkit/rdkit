# $Id$
#
#  Copyright (C) 2016  Novartis Institute of BioMedical Research
#         All Rights Reserved
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
#       products derived from this software without specific prior written permission.
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


from rdkit import Chem
from rdkit.Chem import rdStructChecker

import unittest

data = """310929550
  -OEChem-07211613022D

 55 60  0     0  0  0  0  0  0999 V2000
   -3.6737    1.4194    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8298    2.0868    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7144    1.2798    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5968    2.8875    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8579    0.8673    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    6.4302    1.2798    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.0013   -1.1952    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.5736    0.8673    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.7170    1.2798    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   12.8605    0.8673    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.0093    2.1731    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
   -3.5724   -0.3702    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1434   -0.3702    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.1434    1.2798    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.5724    1.2798    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.8579    0.0423    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2868    0.8673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0013    1.2798    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.2868    0.0423    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7158    0.0423    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.7158    0.8673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0013   -0.3702    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4302   -0.3702    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144   -0.3702    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0423    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4290    0.0423    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.8673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8579    0.0423    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4290    0.8673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144    1.2798    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7144   -1.1952    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7144   -0.3702    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4290    0.8673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.4302   -1.1952    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1447   -0.7827    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.1447    0.0423    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.6077    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7144   -1.1952    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8579    0.8673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4290    0.0423    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.2868    0.8673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1434   -0.3702    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0013    1.2798    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7158    0.8673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.2868    0.0423    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0013   -0.3702    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7158    0.0423    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.1447    0.8673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.2868   -1.6077    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.8592    1.2798    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.2881    1.2798    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.0026    0.8673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.4316    0.8673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.1460    1.2798    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   13.5750    1.2798    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1 11  1  0  0  0  0
  1 17  1  0  0  0  0
  2 11  1  0  0  0  0
  2 18  1  0  0  0  0
  3 27  1  0  0  0  0
  3 33  1  0  0  0  0
  4 11  2  0  0  0  0
  5 28  2  0  0  0  0
  6 44  1  0  0  0  0
  6 48  1  0  0  0  0
  7 46  1  0  0  0  0
  7 49  1  0  0  0  0
  8 50  1  0  0  0  0
  8 51  1  0  0  0  0
  9 52  1  0  0  0  0
  9 53  1  0  0  0  0
 10 54  1  0  0  0  0
 10 55  1  0  0  0  0
 12 19  1  0  0  0  0
 12 28  1  0  0  0  0
 13 26  1  0  0  0  0
 13 28  1  0  0  0  0
 14 33  2  0  0  0  0
 14 39  1  0  0  0  0
 15 39  1  0  0  0  0
 15 41  1  0  0  0  0
 16 39  2  0  0  0  0
 16 42  1  0  0  0  0
 17 18  1  0  0  0  0
 17 19  2  0  0  0  0
 18 21  2  0  0  0  0
 19 22  1  0  0  0  0
 20 21  1  0  0  0  0
 20 22  2  0  0  0  0
 20 23  1  0  0  0  0
 23 34  1  0  0  0  0
 23 35  1  0  0  0  0
 23 36  1  0  0  0  0
 24 25  1  0  0  0  0
 24 26  2  0  0  0  0
 24 31  1  0  0  0  0
 25 27  2  0  0  0  0
 25 32  1  0  0  0  0
 26 29  1  0  0  0  0
 27 30  1  0  0  0  0
 29 30  2  0  0  0  0
 31 37  2  0  0  0  0
 32 38  2  0  0  0  0
 33 40  1  0  0  0  0
 37 38  1  0  0  0  0
 40 42  2  0  0  0  0
 41 43  2  0  0  0  0
 41 45  1  0  0  0  0
 43 44  1  0  0  0  0
 44 47  2  0  0  0  0
 45 46  2  0  0  0  0
 46 47  1  0  0  0  0
 48 50  1  0  0  0  0
 51 52  1  0  0  0  0
 53 54  1  0  0  0  0
M  CHG  1  11   1
M  END
"""


class TestCase(unittest.TestCase):

  def testStructOptions(self):
    rdStructChecker.StructCheckerOptions()

  def testStructChecker(self):
    checker = rdStructChecker.StructChecker()

    m = Chem.MolFromMolBlock(data)
    self.assertTrue(m)

    res = checker.CheckMolStructure(m)
    self.assertEquals(res, rdStructChecker.StructureFlags.ATOM_CHECK_FAILED)


if __name__ == '__main__':
  unittest.main()
