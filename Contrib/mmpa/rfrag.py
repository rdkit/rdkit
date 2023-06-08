# Copyright (c) 2013, GlaxoSmithKline Research & Development Ltd.
# All rights reserved.
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
#     * Neither the name of GlaxoSmithKline Research & Development Ltd.
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
# Created by Jameed Hussain, July 2013
#
# Modifications and optimizations by Greg Landrum, July 2015
#

import re
import sys

from rdkit import Chem
from rdkit.Chem import rdMMPA


def find_correct(f_array):

  core = ""
  side_chains = ""

  for f in f_array:
    attachments = f.count("*")
    if (attachments == 1):
      side_chains = "%s.%s" % (side_chains, f)
    else:
      core = f

  side_chains = side_chains.lstrip('.')

  #cansmi the side chains
  temp = Chem.MolFromSmiles(side_chains)
  side_chains = Chem.MolToSmiles(temp, isomericSmiles=True)

  #and cansmi the core
  temp = Chem.MolFromSmiles(core)
  core = Chem.MolToSmiles(temp, isomericSmiles=True)

  return core, side_chains


def delete_bonds(smi, id, mol, bonds, out):

  #use the same parent mol object and create editable mol
  em = Chem.EditableMol(mol)

  #loop through the bonds to delete
  isotope = 0
  isotope_track = {}
  for i in bonds:
    isotope += 1
    #remove the bond
    em.RemoveBond(i[0], i[1])

    #now add attachment points
    newAtomA = em.AddAtom(Chem.Atom(0))
    em.AddBond(i[0], newAtomA, Chem.BondType.SINGLE)

    newAtomB = em.AddAtom(Chem.Atom(0))
    em.AddBond(i[1], newAtomB, Chem.BondType.SINGLE)

    #keep track of where to put isotopes
    isotope_track[newAtomA] = isotope
    isotope_track[newAtomB] = isotope

  #should be able to get away without sanitising mol
  #as the existing valencies/atoms not changed
  modifiedMol = em.GetMol()

  #canonical smiles can be different with and without the isotopes
  #hence to keep track of duplicates use fragmented_smi_noIsotopes
  fragmented_smi_noIsotopes = Chem.MolToSmiles(modifiedMol, isomericSmiles=True)

  valid = True
  fragments = fragmented_smi_noIsotopes.split(".")

  #check if its a valid triple cut
  if (isotope == 3):
    valid = False
    for f in fragments:
      matchObj = re.search(r'\*.*\*.*\*', f)
      if matchObj:
        valid = True
        break

  if valid:
    if (isotope == 1):
      fragmented_smi_noIsotopes = re.sub(r'\[\*\]', '[*:1]', fragmented_smi_noIsotopes)

      fragments = fragmented_smi_noIsotopes.split(".")

      #print fragmented_smi_noIsotopes
      s1 = Chem.MolFromSmiles(fragments[0])
      s2 = Chem.MolFromSmiles(fragments[1])

      #need to cansmi again as smiles can be different
      output = '%s,%s,,%s.%s' % (smi, id, Chem.MolToSmiles(
        s1, isomericSmiles=True), Chem.MolToSmiles(s2, isomericSmiles=True))
      if output not in out:
        out.add(output)

    elif (isotope >= 2):
      #add the isotope labels
      for key in isotope_track:
        #to add isotope labels
        modifiedMol.GetAtomWithIdx(key).SetIsotope(isotope_track[key])
      fragmented_smi = Chem.MolToSmiles(modifiedMol, isomericSmiles=True)

      #change the isotopes into labels - currently can't add SMARTS or labels to mol
      fragmented_smi = re.sub(r'\[1\*\]', '[*:1]', fragmented_smi)
      fragmented_smi = re.sub(r'\[2\*\]', '[*:2]', fragmented_smi)
      fragmented_smi = re.sub(r'\[3\*\]', '[*:3]', fragmented_smi)

      fragments = fragmented_smi.split(".")

      #identify core/side chains and cansmi them
      core, side_chains = find_correct(fragments)

      #now change the labels on sidechains and core
      #to get the new labels, cansmi the dot-disconnected side chains
      #the first fragment in the side chains has attachment label 1, 2nd: 2, 3rd: 3
      #then change the labels accordingly in the core

      #this is required by the indexing script, as the side-chains are "keys" in the index
      #this ensures the side-chains always have the same numbering

      isotope_track = {}
      side_chain_fragments = side_chains.split(".")

      for s in range(len(side_chain_fragments)):
        matchObj = re.search(r'\[\*\:([123])\]', side_chain_fragments[s])
        if matchObj:
          #add to isotope_track with key: old_isotope, value:
          isotope_track[matchObj.group(1)] = str(s + 1)

      #change the labels if required
      if (isotope_track['1'] != '1'):
        core = re.sub(r'\[\*\:1\]', '[*:XX' + isotope_track['1'] + 'XX]', core)
        side_chains = re.sub(r'\[\*\:1\]', '[*:XX' + isotope_track['1'] + 'XX]', side_chains)
      if (isotope_track['2'] != '2'):
        core = re.sub(r'\[\*\:2\]', '[*:XX' + isotope_track['2'] + 'XX]', core)
        side_chains = re.sub(r'\[\*\:2\]', '[*:XX' + isotope_track['2'] + 'XX]', side_chains)

      if (isotope == 3):
        if (isotope_track['3'] != '3'):
          core = re.sub(r'\[\*\:3\]', '[*:XX' + isotope_track['3'] + 'XX]', core)
          side_chains = re.sub(r'\[\*\:3\]', '[*:XX' + isotope_track['3'] + 'XX]', side_chains)

      #now remove the XX
      core = re.sub('XX', '', core)
      side_chains = re.sub('XX', '', side_chains)

      output = '%s,%s,%s,%s' % (smi, id, core, side_chains)
      if output not in out:
        out.add(output)


def fragment_mol(smi, cid):
  mol = Chem.MolFromSmiles(smi)

  #different cuts can give the same fragments
  #to use outlines to remove them
  outlines = set()

  if mol is None:
    sys.stderr.write("Can't generate mol for: %s\n" % (smi))
  else:
    frags = rdMMPA.FragmentMol(mol, pattern="[#6+0;!$(*=,#[!#6])]!@!=!#[*]", resultsAsMols=False)
    for core, chains in frags:
      output = '%s,%s,%s,%s' % (smi, cid, core, chains)
      if (not (output in outlines)):
        outlines.add(output)
    if not outlines:
      # for molecules with no cuts, output the parent molecule itself
      outlines.add('%s,%s,,' % (smi, cid))

  return outlines


if __name__ == '__main__':

  if (len(sys.argv) >= 2):
    print("Program that fragments a user input set of smiles.")
    print(
      "The program enumerates every single,double and triple acyclic single bond cuts in a molecule.\n"
    )
    print("USAGE: ./rfrag.py <file_of_smiles")
    print("Format of smiles file: SMILES ID (space separated)")
    print("Output: whole mol smiles,ID,core,context\n")
    sys.exit(1)

  #read the STDIN
  for line in sys.stdin:

    line = line.rstrip()

    line_fields = re.split(r'\s|,', line)
    smiles = line_fields[0]
    cmpd_id = line_fields[1]

    #returns a set containing the output
    o = fragment_mol(smiles, cmpd_id)
    for l in o:
      print(l)
"""
optimization work.

Original:
~/RDKit_git/Contrib/mmpa > time head -100 ../../Data/Zinc/zim.smi | python rfrag.py > zim.frags.o

real	0m9.752s
user	0m9.704s
sys	0m0.043s


"""
