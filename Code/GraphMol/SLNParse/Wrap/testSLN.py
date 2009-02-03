#  $Id$
#
#  Copyright (c) 2008, Novartis Institutes for BioMedical Research Inc.
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
#       products derived from this software without specific prior
#       written permission.
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
# Created by Greg Landrum, September 2006
#
from rdkit import Chem
from rdkit.Chem import rdSLNParse
from rdkit import Geometry
from rdkit import RDConfig
import unittest
import os,sys

class TestCase(unittest.TestCase) :
  def setUp(self):
    self.dataDir = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','SLNParse','testData')

  def test1Basics(self):
    m1 = rdSLNParse.MolFromSLN('CH3CH3')
    self.failUnless(m1)
    self.failUnless(m1.GetNumAtoms()==2)

    m1 = rdSLNParse.MolFromSLN('C[1]H:CH:CH:CH:CH:CH:@1')
    self.failUnless(m1)
    self.failUnless(m1.GetNumAtoms()==6)
    
  def test2Queries(self):
    patt = rdSLNParse.MolFromQuerySLN('C[HC=2]~O')
    self.failUnless(patt)
    self.failUnless(patt.GetNumAtoms()==2)
    
    m=Chem.MolFromSmiles('COCC=O')
    self.failUnless(m.HasSubstructMatch(patt))
    ms = m.GetSubstructMatches(patt)
    self.failUnless(len(ms)==1)
    
if __name__ == '__main__':
  unittest.main()
