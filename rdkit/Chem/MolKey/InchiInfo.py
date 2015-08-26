#
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
# Retrieve stereo and tautomer information from the the InChI string
# Created on Sep 23, 2010
# Original author: Thomas Muellerk muelleth
import logging
import re
import unittest
from rdkit import Chem

from rdkit.Chem import inchi
if not inchi.INCHI_AVAILABLE:
    raise ImportError("This code requires the RDKit to be built with InChI suport")

def _is_achiral_by_symmetry(INCHI) :
    mol = Chem.MolFromInchi(INCHI)
    if not mol :
        mol = Chem.MolFromInchi('InChI=1/{0}'.format(INCHI))
        
    try :
        list_chiral = Chem.FindMolChiralCenters(mol, True, True)
    except Exception :
        return False
 
#   is there any real chiral centre?
    return len(list_chiral) == 0

console = logging.StreamHandler()
UPD_APP = logging.getLogger('inchiinfo.application')   # application runtime information 

version_re = re.compile('(.*?)/(.*)')  # get version
reconnected_re = re.compile('(.*?)/r(.*)')  # reconnected layer?
fixed_h_re = re.compile('(.*?)/f(.*)')  # fixed-H layer?
isotope_re = re.compile('(.*?)/i(.*)')  # isotope layer?

stereo_re = re.compile('.*\/t(.*?)\/.*')
stereo_all_re = re.compile('.*\/t([^\/]+)')
undef_stereo_re = re.compile('(\d+)\?')
all_stereo_re = re.compile('(\d+)[?+-]')
defined_stereo_re = re.compile('(\d+)[+-]')
h_layer_re = re.compile('.*\/h(.*)\/?')
mobile_h_group_re = re.compile('(\(H.+?\))')
mobile_h_atoms_re = re.compile(',(\d+)')
 
class InchiInfo(object): 
    
    def __init__(self, inchi_str):
        (version, rest) = version_re.match(inchi_str).groups()
        reconn_match = reconnected_re.match(rest)
                    
        connection_layers = {}
        if reconn_match:
            (connection_layers['id_disconnected'], connection_layers['id_reconnected']) = reconn_match.groups()
        else:
            (connection_layers['id']) = rest
        
        fixed_h_layers = {}
        for conn_layer in connection_layers:
            fixed_h_layers[conn_layer] = {}
            fixed_match = fixed_h_re.match(connection_layers[conn_layer])
            if fixed_match:
                (fixed_h_layers[conn_layer]['main'], fixed_h_layers[conn_layer]['fixed_h']) = fixed_match.groups()
            else:    
                fixed_h_layers[conn_layer]['main'] = connection_layers[conn_layer]

        inchi = {}
        for i0_layer in fixed_h_layers:
            inchi[i0_layer] = {}
            for i1_layer in fixed_h_layers[i0_layer]:
                inchi[i0_layer][i1_layer] = {}
                iso_match = isotope_re.match(fixed_h_layers[i0_layer][i1_layer])
                if iso_match:
                    (inchi[i0_layer][i1_layer]['non-isotopic'], inchi[i0_layer][i1_layer]['isotopic']) = iso_match.groups() 
                else:
                    inchi[i0_layer][i1_layer]['non-isotopic'] = fixed_h_layers[i0_layer][i1_layer]
            
        self.parsed_inchi = inchi

    def get_sp3_stereo(self):
        ''' retrieve sp3 stereo information
        return a 4-item tuple containing
        1) Number of stereocenters detected. If 0, the remaining items of the tuple = None
        2) Number of undefined stereocenters. Must be smaller or equal to above
        3) True if the molecule is a meso form (with chiral centers and a plane of symmetry)
        4) Comma-separated list of internal atom numbers with sp3 stereochemistry
        ''' 
        sp3_stereo = {}

        for con_layer in self.parsed_inchi:
            for fixed_layer in self.parsed_inchi[con_layer]:
                sp3_stereo[fixed_layer] = {}
                for iso_layer in self.parsed_inchi[con_layer][fixed_layer]:
                    sp3_stereo[fixed_layer][iso_layer] = {}
                    stereo_match = stereo_re.match(self.parsed_inchi[con_layer][fixed_layer][iso_layer])
                    stereo_all_match = stereo_all_re.match(self.parsed_inchi[con_layer][fixed_layer][iso_layer])
                    num_stereo = 0
                    num_undef_stereo = 0
                    is_meso = False
                    stereo = ''
                    stereo_centers = []
                    undef_stereo_centers = []
                    #    match patterns with defined and undefined stereo
                    if stereo_match:
                        stereo = stereo_match.group(1)
                    #    match patterns with only undefined stereo or for the MESO case
                    elif stereo_all_match :
                        stereo = stereo_all_match.group(1)
                        is_meso = len(defined_stereo_re.findall(stereo)) > 1
                    #    Number of ALL stereo centres
                    stereo_centers = all_stereo_re.findall(stereo)
                    num_stereo = len(stereo_centers)
                    undef_stereo_centers = undef_stereo_re.findall(stereo)
                    num_undef_stereo = len(undef_stereo_centers)
                    #    Meso centres    --    VT     --    2011.12.08
                    inchi_layer = self.parsed_inchi[con_layer][fixed_layer][iso_layer]
                    is_meso = is_meso or (num_undef_stereo > 1 and _is_achiral_by_symmetry(inchi_layer))
                    sp3_stereo[fixed_layer][iso_layer] = (num_stereo, num_undef_stereo, is_meso, stereo)
        return sp3_stereo                    

    def get_mobile_h(self):
        ''' retrieve mobile H (tautomer) information
        return a 2-item tuple containing
        1) Number of mobile hydrogen groups detected. If 0, next item = '' 
        2) List of groups   
        '''
        mobile_h = {}
        for con_layer in self.parsed_inchi:
            for fixed_layer in self.parsed_inchi[con_layer]:
                mobile_h[fixed_layer] = {}
                for iso_layer in self.parsed_inchi[con_layer][fixed_layer]:
                    num_groups = 0
                    mobile_h_groups = ''    
                    h_layer_match = h_layer_re.match(self.parsed_inchi[con_layer][fixed_layer][iso_layer])
                    if h_layer_match:
                        mobile_h_matches = mobile_h_group_re.findall(h_layer_match.group(1))
                        num_groups = len(mobile_h_matches)
                        mobile_h_groups = ','.join(mobile_h_matches)
                    mobile_h[fixed_layer][iso_layer] = (num_groups, mobile_h_groups)
        return mobile_h            

# Test molecules as InChI strings
# tautomers
GUANINE='InChI=1S/C5H5N5O/c6-5-9-3-2(4(11)10-5)7-1-8-3/h1H0,(H4,6,7,8,9,10,11)'
# 'N=C(-O)N', '/FixedH /SUU'
UREA1 = 'InChI=1/CH4N2O/c2-1(3)4/h(H4,2,3,4)/f/h2,4H,3H2/b2-1?'
# 'NC(=O)N', '/FixedH /SUU'
UREA2 = 'InChI=1/CH4N2O/c2-1(3)4/h(H4,2,3,4)/f/h2-3H2'
TRITIATED_UREA='InChI=1S/CH4N2O/c2-1(3)4/h(H4,2,3,4)/i/hT3'
DEUTERATED_UREA='InChI=1S/CH4N2O/c2-1(3)4/h(H4,2,3,4)/i/hD2'
ACETIC_ACID='InChI=1S/C3H6O2/c1-2-3(4)5/h2H2,1H3,(H,4,5)'
ACETATE='InChI=1S/C3H6O2/c1-2-3(4)5/h2H2,1H3,(H,4,5)/p-1'
mobile1='InChI=1S/C5H5N3O2/c6-4(9)3-1-7-2-8-5(3)10/h1-2H,(H2,6,9)(H,7,8,10)' # invented
mobile2='InChI=1S/C7H10N4O/c1-4-2-5(3-6(8)12)11-7(9)10-4/h2H,3H2,1H3,(H2,8,12)(H2,9,10,11)'

# sp3 stereo
sugar1='InChI=1S/C14H20O9/c1-6-11(20-7(2)15)12(21-8(3)16)13(22-9(4)17)14(19-6)23-10(5)18/h6,11-14H,1-5H3/t6-,11-,12+,13+,14?/m0/s1' # L-rhamnopyranose (source: chemspider)
sugar2='InChI=1S/C12H20O6/c1-11(2)14-5-6(16-11)8-7(13)9-10(15-8)18-12(3,4)17-9/h6-10,13H,5H2,1-4H3/t6-,7-,8-,9-,10-/m1/s1' # MFCD00135634 (Diacetone-D-Glucose, souce: chemspider)
sp3_unk='InChI=1S/C12H21NO4/c1-8(2)10(12(15)16-3)13-11(14)9-5-4-6-17-7-9/h8-10H,4-7H2,1-3H3,(H,13,14)/t9?,10-/m0/s1' # derived from ChemSpider 34044335

class TestInchiInfo(unittest.TestCase):

    def doTest(self, inchi, numSp3=0, numUndefSp3=0, numMobileHGroups=0, layer='non-isotopic'):
        ii = InchiInfo(inchi)
        (nSp3, nUndefSp3, isMeso, sp3Atoms) = ii.get_sp3_stereo()['main'][layer]
        self.assertEqual(nSp3, numSp3)
        self.assertEqual(nUndefSp3, numUndefSp3)
            
        (nMobileHGroups, mobileHGroups) = ii.get_mobile_h()['main'][layer]    
        self.assertEqual(nMobileHGroups, numMobileHGroups)    
        
    def testGuanine(self): 
      self.doTest(GUANINE, 0, 0, 1) 
    def testTritiatedUrea(self): 
      self.doTest(TRITIATED_UREA, 0, 0, 1)
    def testDeuteratedUrea(self): 
      self.doTest(DEUTERATED_UREA, 0, 0, 1)
    def testAceticAcid(self): 
      self.doTest(ACETIC_ACID, 0, 0, 1)
    def testAcetate(self): 
      self.doTest(ACETATE, 0, 0, 1)
    
    def testMobile1(self): 
      self.doTest(mobile1, 0, 0, 2)
    def testMobile2(self): 
      self.doTest(mobile2, 0, 0, 2)

    
    # sp3 stereo
    def testSugar1(self): 
      self.doTest(sugar1, 5, 1, 0)
    def testSugar2(self): 
      self.doTest(sugar2, 5, 0, 0)
    def testSP3_unk(self): 
      self.doTest(sp3_unk, 2, 1, 1)

if __name__ == '__main__':
    unittest.main()
    
    
