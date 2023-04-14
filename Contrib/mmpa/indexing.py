# Copyright (c) 2012, GlaxoSmithKline Research & Development Ltd.
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
# Created by Jameed Hussain, September 2012

import re
import sys
from optparse import OptionParser

from rdkit import Chem


def heavy_atom_count(smi):

  m = Chem.MolFromSmiles(smi)
  return m.GetNumAtoms()


def add_to_index(smi, attachments, cmpd_heavy):

  result = False
  core_size = heavy_atom_count(smi) - attachments

  if (use_ratio):
    core_ratio = float(core_size) / float(cmpd_heavy)
    if (core_ratio <= ratio):
      result = True
  else:
    if (core_size <= max_size):
      result = True

  return result


def get_symmetry_class(smi):
  symmetry = []

  m = Chem.MolFromSmiles(smi)
  symmetry_classes = Chem.CanonicalRankAtoms(m, breakTies=False)

  #get the symmetry class of the attachments points
  #Note: 1st star is the zero index,
  #2nd star is first index, etc
  for atom, symmetry_class in zip(m.GetAtoms(), symmetry_classes):
    if (atom.GetMass() == 0):
      symmetry.append(symmetry_class)

  return symmetry


def cansmirk(lhs, rhs, context):

  #cansmirk algorithm
  #1) cansmi the LHS.
  #2) For the LHS the 1st star will have label 1, 2nd star will have label 2 and so on
  #3) Do a symmetry check of lhs and rhs and use that to decide if the labels on
  #   RHS or/and context need to change.
  #4) For the rhs, if you have a choice (ie. two attachment points are symmetrically
  #   equivalent), always put the label with lower numerical value on the earlier
  #   attachment point on the cansmi-ed smiles

  #print "in: %s,%s" % (lhs,rhs)

  isotope_track = {}
  #if the star count of lhs/context/rhs is 1, single cut
  stars = lhs.count("*")

  if (stars > 1):
    #get the symmetry class of stars of lhs and rhs
    lhs_sym = get_symmetry_class(lhs)
    rhs_sym = get_symmetry_class(rhs)

  #deal with double cuts
  if (stars == 2):
    #simple cases
    #unsymmetric lhs and unsymmetric rhs
    if ((lhs_sym[0] != lhs_sym[1]) and (rhs_sym[0] != rhs_sym[1])):
      #get 1st and 2nd labels and store the new label for it in isotope_track
      #structure: isotope_track[old_label]=new_label (as strings)
      isotope_track = build_track_dictionary(lhs, stars)

      #switch labels using isotope track
      lhs = switch_labels_on_position(lhs)
      rhs = switch_labels(isotope_track, stars, rhs)
      context = switch_labels(isotope_track, stars, context)

    #symmetric lhs and symmetric rhs
    elif ((lhs_sym[0] == lhs_sym[1]) and (rhs_sym[0] == rhs_sym[1])):
      #the points are all equivalent so change labels on lhs and rhs based on position
      #labels on context don't need to change
      lhs = switch_labels_on_position(lhs)
      rhs = switch_labels_on_position(rhs)

    #more difficult cases..
    #symmetric lhs and unsymmetric rhs
    elif ((lhs_sym[0] == lhs_sym[1]) and (rhs_sym[0] != rhs_sym[1])):
      #switch labels lhs based on position
      lhs = switch_labels_on_position(lhs)
      #change labels on rhs based on position but need to record
      #the changes as need to apply them to the context
      isotope_track = build_track_dictionary(rhs, stars)
      rhs = switch_labels_on_position(rhs)
      context = switch_labels(isotope_track, stars, context)

    #unsymmetric lhs and symmetric rhs
    elif ((lhs_sym[0] != lhs_sym[1]) and (rhs_sym[0] == rhs_sym[1])):
      #change labels on lhs based on position but need to record
      #the changes as need to apply them to the context
      isotope_track = build_track_dictionary(lhs, stars)
      lhs = switch_labels_on_position(lhs)
      context = switch_labels(isotope_track, stars, context)
      #as rhs is symmetric, positions are equivalent so change labels on position
      rhs = switch_labels_on_position(rhs)

    #deal with triple cut
    #unwieldy code but most readable I can make it
  elif (stars == 3):
    #simple cases
    #completely symmetric lhs and completely symmetric rhs
    if (((lhs_sym[0] == lhs_sym[1]) and (lhs_sym[1] == lhs_sym[2]) and (lhs_sym[0] == lhs_sym[2]))
        and ((rhs_sym[0] == rhs_sym[1]) and (rhs_sym[1] == rhs_sym[2]) and
             (rhs_sym[0] == rhs_sym[2]))):
      #the points are all equivalent so change labels on lhs and rhs based on position
      #labels on context don't need to change
      lhs = switch_labels_on_position(lhs)
      rhs = switch_labels_on_position(rhs)

    #completely symmetric lhs and completely unsymmetric rhs
    elif (((lhs_sym[0] == lhs_sym[1]) and (lhs_sym[1] == lhs_sym[2]) and (lhs_sym[0] == lhs_sym[2]))
          and ((rhs_sym[0] != rhs_sym[1]) and (rhs_sym[1] != rhs_sym[2]) and
               (rhs_sym[0] != rhs_sym[2]))):

      #alter lhs in usual way
      lhs = switch_labels_on_position(lhs)
      #change labels on rhs based on position but need to record
      #the changes as need to apply them to the context
      isotope_track = build_track_dictionary(rhs, stars)
      rhs = switch_labels_on_position(rhs)
      context = switch_labels(isotope_track, stars, context)

    #completely unsymmetric lhs and completely unsymmetric rhs
    elif (((lhs_sym[0] != lhs_sym[1]) and (lhs_sym[1] != lhs_sym[2]) and (lhs_sym[0] != lhs_sym[2]))
          and ((rhs_sym[0] != rhs_sym[1]) and (rhs_sym[1] != rhs_sym[2]) and
               (rhs_sym[0] != rhs_sym[2]))):

      #build the isotope track
      isotope_track = build_track_dictionary(lhs, stars)
      #alter lhs in usual way
      lhs = switch_labels_on_position(lhs)
      #change rhs and context based on isotope_track
      rhs = switch_labels(isotope_track, stars, rhs)
      context = switch_labels(isotope_track, stars, context)

    #completely unsymmetric lhs and completely symmetric rhs
    elif (((lhs_sym[0] != lhs_sym[1]) and (lhs_sym[1] != lhs_sym[2]) and (lhs_sym[0] != lhs_sym[2]))
          and ((rhs_sym[0] == rhs_sym[1]) and (rhs_sym[1] == rhs_sym[2]) and
               (rhs_sym[0] == rhs_sym[2]))):

      #build isotope trach on lhs
      isotope_track = build_track_dictionary(lhs, stars)
      #alter lhs in usual way
      lhs = switch_labels_on_position(lhs)
      #change labels on context
      context = switch_labels(isotope_track, stars, context)
      #all positions on rhs equivalent so add labels on position
      rhs = switch_labels_on_position(rhs)

    #more difficult cases, partial symmetry
    #completely unsymmetric on lhs and partial symmetry on rhs
    elif ((lhs_sym[0] != lhs_sym[1]) and (lhs_sym[1] != lhs_sym[2]) and (lhs_sym[0] != lhs_sym[2])):
      #build the isotope track
      isotope_track = build_track_dictionary(lhs, stars)
      #alter lhs in usual way
      lhs = switch_labels_on_position(lhs)
      #change rhs and context based on isotope_track
      rhs = switch_labels(isotope_track, stars, rhs)
      context = switch_labels(isotope_track, stars, context)

      #tweak positions on rhs based on symmetry
      #rhs 1,2 equivalent
      if (rhs_sym[0] == rhs_sym[1]):
        #tweak rhs position 1 and 2 as they are symmetric
        rhs = switch_specific_labels_on_symmetry(rhs, rhs_sym, 1, 2)

      #rhs 2,3 equivalent
      elif (rhs_sym[1] == rhs_sym[2]):
        #tweak rhs position 1 and 2 as they are symmetric
        rhs = switch_specific_labels_on_symmetry(rhs, rhs_sym, 2, 3)

      #rhs 1,3 equivalent - try for larger set in future
      elif (rhs_sym[0] == rhs_sym[2]):
        #tweak rhs position 1 and 2 as they are symmetric
        rhs = switch_specific_labels_on_symmetry(rhs, rhs_sym, 1, 3)

      #now we are left with things with partial symmetry on lhs and not completely symmetric or unsymmetric on rhs
    else:
      #lhs 1,2,3 equivalent and any sort of partial symmetry on rhs
      if ((lhs_sym[0] == lhs_sym[1]) and (lhs_sym[1] == lhs_sym[2]) and (lhs_sym[0] == lhs_sym[2])):

        #alter lhs in usual way
        lhs = switch_labels_on_position(lhs)
        #change labels on rhs based on position but need to record
        #the changes as need to apply them to the context
        isotope_track = build_track_dictionary(rhs, stars)
        rhs = switch_labels_on_position(rhs)
        context = switch_labels(isotope_track, stars, context)

      #now deal partial symmetry on lhs or rhs.
      #Cases where:
      #lhs 1,2 equivalent
      #lhs 2,3 equivalent
      #lhs 1,3 equivalent
      else:
        #build isotope track on lhs
        isotope_track = build_track_dictionary(lhs, stars)
        #alter lhs in usual way
        lhs = switch_labels_on_position(lhs)
        #change rhs and context based on isotope_track
        rhs = switch_labels(isotope_track, stars, rhs)
        context = switch_labels(isotope_track, stars, context)

        #tweak positions on rhs based on symmetry

        #lhs 1,2 equivalent
        if (lhs_sym[0] == lhs_sym[1]):
          #tweak rhs position 1 and 2 as they are symmetric on lhs
          rhs = switch_specific_labels_on_symmetry(rhs, rhs_sym, 1, 2)

        #lhs 2,3 equivalent
        elif (lhs_sym[1] == lhs_sym[2]):
          #tweak rhs position 1 and 2 as they are symmetric on lhs
          rhs = switch_specific_labels_on_symmetry(rhs, rhs_sym, 2, 3)

        #lhs 1,3 equivalent - try for larger set in future
        elif (lhs_sym[0] == lhs_sym[2]):
          #tweak rhs position 1 and 2 as they are symmetric on lhs
          rhs = switch_specific_labels_on_symmetry(rhs, rhs_sym, 1, 3)

  smirk = "%s>>%s" % (lhs, rhs)

  return smirk, context


def switch_specific_labels_on_symmetry(smi, symmetry_class, a, b):

  #check if a and b positions are symmetrically equivalent
  #if equivalent, swap labels if the lower numerical label is not on the
  #1st symmetrically equivalent attachment points in the smi

  if (symmetry_class[a - 1] == symmetry_class[b - 1]):
    #what are the labels on a and b

    matchObj = re.search(r'\[\*\:([123])\].*\[\*\:([123])\].*\[\*\:([123])\]', smi)
    if matchObj:
      #if the higher label comes first, fix
      if (int(matchObj.group(a)) > int(matchObj.group(b))):
        #if(int(matchObj.group(1)) > int(matchObj.group(2))):
        smi = re.sub(r'\[\*\:' + matchObj.group(a) + r'\]', '[*:XX' + matchObj.group(b) + 'XX]',
                     smi)
        smi = re.sub(r'\[\*\:' + matchObj.group(b) + r'\]', '[*:XX' + matchObj.group(a) + 'XX]',
                     smi)
        smi = re.sub('XX', '', smi)

  return smi


def switch_labels_on_position(smi):

  #move the labels in order of position
  smi = re.sub(r'\[\*\:[123]\]', '[*:XX1XX]', smi, 1)
  smi = re.sub(r'\[\*\:[123]\]', '[*:XX2XX]', smi, 1)
  smi = re.sub(r'\[\*\:[123]\]', '[*:XX3XX]', smi, 1)
  smi = re.sub('XX', '', smi)

  return smi


def switch_labels(track, stars, smi):

  #switch labels based on the input dictionary track
  if (stars > 1):
    #for k in track:
    #        print "old: %s, new: %s" % (k,track[k])

    if (track['1'] != '1'):
      smi = re.sub(r'\[\*\:1\]', '[*:XX' + track['1'] + 'XX]', smi)

    if (track['2'] != '2'):
      smi = re.sub(r'\[\*\:2\]', '[*:XX' + track['2'] + 'XX]', smi)

    if (stars == 3):
      if (track['3'] != '3'):
        smi = re.sub(r'\[\*\:3\]', '[*:XX' + track['3'] + 'XX]', smi)

    #now remove the XX
    smi = re.sub('XX', '', smi)

  return smi


def build_track_dictionary(smi, stars):

  isotope_track = {}

  #find 1st label, record it in isotope_track as key, with value being the
  #new label based on its position (1st star is 1, 2nd star 2 etc.)
  if (stars == 2):
    matchObj = re.search(r'\[\*\:([123])\].*\[\*\:([123])\]', smi)
    if matchObj:
      isotope_track[matchObj.group(1)] = '1'
      isotope_track[matchObj.group(2)] = '2'

  elif (stars == 3):
    matchObj = re.search(r'\[\*\:([123])\].*\[\*\:([123])\].*\[\*\:([123])\]', smi)
    if matchObj:
      isotope_track[matchObj.group(1)] = '1'
      isotope_track[matchObj.group(2)] = '2'
      isotope_track[matchObj.group(3)] = '3'

  return isotope_track


def index_hydrogen_change():
  #Algorithm details
  #have an index of common fragment(key) => fragments connected to it (values)
  #Need to add *-H to the values where appropriate - and its
  #appropriate when the key is what you would get if you chopped a H off a cmpd.
  #Therefore simply need to check if key with the * replaced with a H is
  #the same as any full smiles in the set
  #
  #Specific details:
  #1) Loop through keys of index
  #2) If key is the result of a single cut (so contains only 1 *) replace the * with H, and cansmi
  #3) If full smiles matches key in hash above, add *-H to that fragment index.

  for key in index:

    attachments = key.count('*')
    #print attachments

    if (attachments == 1):

      smi = key

      #simple method
      smi = re.sub(r'\[\*\:1\]', '[H]', smi)

      #now cansmi it
      temp = Chem.MolFromSmiles(smi)

      if temp is None:
        sys.stderr.write('Error with key: %s, Added H: %s\n' % (key, smi))
      else:
        c_smi = Chem.MolToSmiles(temp, isomericSmiles=True)

        if (c_smi in smi_to_id):
          core = "[*:1][H]"
          id = smi_to_id[c_smi]

          value = "%s;t%s" % (id, core)
          #add to index
          index[key].append(value)


if __name__ == '__main__':

  #note max heavy atom count does not
  #include the attachment points (*)
  max_size = 10
  ratio = 0.3
  use_ratio = False

  index = {}
  smi_to_id = {}
  id_to_smi = {}

  id_to_heavy = {}

  #set up the command line options
  #parser = OptionParser()
  parser = OptionParser(description="Program to generate MMPs")
  parser.add_option(
    '-s', '--symmetric', default=False, action='store_true', dest='sym', help=
    'Output symmetrically equivalent MMPs, i.e output both cmpd1,cmpd2, SMIRKS:A>>B and cmpd2,cmpd1, SMIRKS:B>>A'
  )
  parser.add_option(
    '-m', '--maxsize', action='store', dest='maxsize', type='int', help=
    'Maximum size of change (in heavy atoms) allowed in matched molecular pairs identified. DEFAULT=10. \
                      Note: This option overrides the ratio option if both are specified.')
  parser.add_option(
    '-r', '--ratio', action='store', dest='ratio', type='float', help=
    'Maximum ratio of change allowed in matched molecular pairs identified. The ratio is: size of change / \
                      size of cmpd (in terms of heavy atoms). DEFAULT=0.3. Note: If this option is used with the maxsize option, the maxsize option will be used.'
  )

  #parse the command line options
  (options, args) = parser.parse_args()

  #print options
  if options.maxsize is not None:
    max_size = options.maxsize
  elif options.ratio is not None:
    ratio = options.ratio
    if (ratio >= 1):
      print("Ratio specified: %s. Ratio needs to be less than 1.")
      sys.exit(1)
    use_ratio = True

  #read the STDIN
  for line in sys.stdin:

    line = line.rstrip()
    smi, id, core, context = line.split(',')

    #fill in dictionaries
    smi_to_id[smi] = id
    id_to_smi[id] = smi

    #if using the ratio option, check if heavy atom
    #of mol already calculated. If not, calculate and store
    cmpd_heavy = None
    if (use_ratio):
      if id not in id_to_heavy:
        id_to_heavy[id] = heavy_atom_count(smi)

      cmpd_heavy = id_to_heavy[id]

    #deal with cmpds that have not been fragmented
    if (len(core) == 0) and (len(context) == 0):
      continue

    #deal with single cuts
    if (len(core) == 0):
      side_chains = context.split('.')

      #minus 1 for the attachment pt
      if add_to_index(side_chains[1], 1, cmpd_heavy):
        context = side_chains[0]
        core = side_chains[1]

        value = "%s;t%s" % (id, core)

        #add the array if no key exists
        #add the context with id to index
        index.setdefault(context, []).append(value)

      #minus 1 for the attachment pt
      if add_to_index(side_chains[0], 1, cmpd_heavy):
        context = side_chains[1]
        core = side_chains[0]

        value = "%s;t%s" % (id, core)

        #add the array if no key exists
        #add the context with id to index
        index.setdefault(context, []).append(value)

    #double or triple cut
    else:

      attachments = core.count('*')

      if add_to_index(core, attachments, cmpd_heavy):
        value = "%s;t%s" % (id, core)

        #add the array if no key exists
        #add the context with id to index
        index.setdefault(context, []).append(value)

  #index the H change
  index_hydrogen_change()

  #Now index is ready

  #loop through the index
  for key in index:

    total = len(index[key])

    #check if have more than one value
    if (total == 1):
      continue

    for xa in range(total):

      for xb in range(xa, total):

        if (xa != xb):
          #now generate the pairs

          id_a, core_a = index[key][xa].split(";t")
          id_b, core_b = index[key][xb].split(";t")

          #make sure pairs are not same molecule
          if (id_a != id_b):

            #make sure LHS and RHS of SMIRKS are not the same
            if (core_a != core_b):

              smirks, context = cansmirk(core_a, core_b, key)
              print("%s,%s,%s,%s,%s,%s" %
                    (id_to_smi[id_a], id_to_smi[id_b], id_a, id_b, smirks, context))

              #deal with symmetry switch
              if options.sym:
                smirks, context = cansmirk(core_b, core_a, key)
                print("%s,%s,%s,%s,%s,%s" %
                      (id_to_smi[id_b], id_to_smi[id_a], id_b, id_a, smirks, context))
