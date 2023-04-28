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
# Retrieve stereo and tautomer information from the InChI string
# Created on Sep 23, 2010
# Original author: Thomas Muellerk muelleth
import logging
import re

from rdkit import Chem
from rdkit.Chem import inchi

if not inchi.INCHI_AVAILABLE:
  raise ImportError("This code requires the RDKit to be built with InChI support")


def _is_achiral_by_symmetry(INCHI):
  mol = Chem.MolFromInchi(INCHI)
  if not mol:
    mol = Chem.MolFromInchi(f'InChI=1/{INCHI}')
  try:
    # is there any real chiral centre?
    return len(Chem.FindMolChiralCenters(mol, True, True)) == 0
  except Exception:
    return False


console = logging.StreamHandler()
UPD_APP = logging.getLogger('inchiinfo.application')  # application runtime information

version_re = re.compile('(.*?)/(.*)')  # get version
reconnected_re = re.compile('(.*?)/r(.*)')  # reconnected layer?
fixed_h_re = re.compile('(.*?)/f(.*)')  # fixed-H layer?
isotope_re = re.compile('(.*?)/i(.*)')  # isotope layer?

stereo_re = re.compile(r'.*/t(.*?)/.*')
stereo_all_re = re.compile(r'.*/t([^/]+)')
undef_stereo_re = re.compile(r'(\d+)\?')
all_stereo_re = re.compile(r'(\d+)[?+-]')
defined_stereo_re = re.compile(r'(\d+)[+-]')
h_layer_re = re.compile(r'.*/h(.*)/?')
mobile_h_group_re = re.compile(r'(\(H.+?\))')
mobile_h_atoms_re = re.compile(r',(\d+)')


class InchiInfo(object):

  def __init__(self, inchi_str):
    _, rest = version_re.match(inchi_str).groups()
    reconn_match = reconnected_re.match(rest)

    connection_layers = {}
    if reconn_match:
      connection_layers['id_disconnected'], connection_layers[
        'id_reconnected'] = reconn_match.groups()
    else:
      connection_layers['id'] = rest

    fixed_h_layers = {}
    for conn_layer in connection_layers:
      fixed_h_layers[conn_layer] = {}
      fixed_match = fixed_h_re.match(connection_layers[conn_layer])
      if fixed_match:
        fixed_h_layers[conn_layer]['main'], fixed_h_layers[conn_layer][
          'fixed_h'] = fixed_match.groups()
      else:
        fixed_h_layers[conn_layer]['main'] = connection_layers[conn_layer]

    inchi = {}
    for i0_layer in fixed_h_layers:
      inchi[i0_layer] = {}
      for i1_layer in fixed_h_layers[i0_layer]:
        inchi[i0_layer][i1_layer] = {}
        iso_match = isotope_re.match(fixed_h_layers[i0_layer][i1_layer])
        if iso_match:
          inchi[i0_layer][i1_layer]['non-isotopic'], inchi[i0_layer][i1_layer][
            'isotopic'] = iso_match.groups()
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
          stereo_all_match = stereo_all_re.match(
            self.parsed_inchi[con_layer][fixed_layer][iso_layer])

          num_stereo = 0
          num_undef_stereo = 0
          is_meso = False
          stereo = ''
          #    match patterns with defined and undefined stereo
          if stereo_match:
            stereo = stereo_match.group(1)
          #    match patterns with only undefined stereo or for the MESO case
          elif stereo_all_match:
            stereo = stereo_all_match.group(1)
            is_meso = len(defined_stereo_re.findall(stereo)) > 1

          #    Number of ALL stereo centres
          num_stereo = len(all_stereo_re.findall(stereo))
          num_undef_stereo = len(undef_stereo_re.findall(stereo))

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
