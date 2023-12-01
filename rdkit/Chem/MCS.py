# This work was funded by Roche and generously donated to the free
# and open source cheminformatics community.
import warnings

warnings.simplefilter('default', DeprecationWarning)
warnings.warn("The rdkit.Chem.MCS module is deprecated; please use rdkit.Chem.rdFMCS instead.",
              DeprecationWarning, stacklevel=2)
## Copyright (c) 2012 Andrew Dalke Scientific AB
## Andrew Dalke <dalke@dalkescientific.com>
##
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##
##   * Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##
##   * Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
from rdkit.Chem import fmcs
from rdkit.Chem.fmcs import Default

"""MCS - find a Maximum Common Substructure

This software finds the maximum common substructure of a set of
structures and reports it as a SMARTS string.

The SMARTS string depends on the desired match properties. For
example, if ring atoms are only allowed to match ring atoms then an
aliphatic ring carbon in the query is converted to the SMARTS "[C;R]",
and the double-bond ring bond converted to "=;@" while the respective
chain-only version are "[C;!R]" and "=;!@".

"""

# The simplified algorithm description is:
#
#   best_substructure = None
#   pick one structure as the query, and other as the targets
#   for each substructure in the query graph:
#     convert it to a SMARTS string based on the desired match properties
#     if the SMARTS pattern exists in all of the targets:
#        then this is a common substructure
#        keep track of the maximum such common structure,
#
# The algorithm will usually take a long time. There are several
# ways to speed it up.
#
# == Bond elimination ==
#
# As the first step, remove bonds which obviously cannot be part of the
# MCS.
#
# This requires atom and bond type information, which I store as SMARTS
# patterns. A bond can only be in the MCS if its canonical bond type is
# present in all of the structures. A bond type is string made of the
# SMARTS for one atom, the SMARTS for the bond, and the SMARTS for the
# other atom. The canonical bond type is the lexographically smaller of
# the two possible bond types for a bond.
#
# The atom and bond SMARTS depend on the type comparison used.
#
# The "ring-matches-ring-only" option adds an "@" or "!@" to the bond
# SMARTS, so that the canonical bondtype for "C-C" becomes [#6]-@[#6] or
# [#6]-!@[#6] if the bond is in a ring or not in a ring, and if atoms
# are compared by element and bonds are compared by bondtype. (This
# option does not add "R" or "!R" to the atom SMARTS because there
# should be a single bond in the MCS of c1ccccc1O and CO.)
#
# The result of all of this atom and bond typing is a "TypedMolecule"
# for each input structure.
#
# I then find which canonical bondtypes are present in all of the
# structures. I convert each TypedMolecule into a
# FragmentedTypedMolecule which has the same atom information but only
# those bonds whose bondtypes are in all of the structures. This can
# break a structure into multiple, disconnected fragments, hence the
# name.
#
# (BTW, I would like to use the fragmented molecules as the targets
# because I think the SMARTS match would go faster, but the RDKit SMARTS
# matcher doesn't like them. I think it's because the new molecule
# hasn't been sanitized and the underlying data structure the ring
# information doesn't exist. Instead, I use the input structures for the
# SMARTS match.)
#
# == Use the structure with the smallest largest fragment as the query ==
# == and sort the targets by the smallest largest fragment             ==
#
# I pick one of the FragmentedTypedMolecule instances as the source of
# substructure enumeration. Which one?
#
# My heuristic is to use the one with the smallest largest fragment.
# Hopefully it produces the least number of subgraphs, but that's also
# related to the number of rings, so a large linear graph will product
# fewer subgraphs than a small fused ring system. I don't know how to
# quantify that.
#
# For each of the fragmented structures, I find the number of atoms in
# the fragment with the most atoms, and I find the number of bonds in
# the fragment with the most bonds. These might not be the same
# fragment.
#
# I sort the input structures by the number of bonds in the largest
# fragment, with ties broken first on the number of atoms, and then on
# the input order. The smallest such structure is the query structure,
# and the remaining are the targets.
#
# == Use a breadth-first search and a priority queue to    ==
# == enumerate the fragment subgraphs                      ==
#
# I extract each of the fragments from the FragmentedTypedMolecule into
# a TypedFragment, which I use to make an EnumerationMolecule. An
# enumeration molecule contains a pair of directed edges for each atom,
# which simplifies the enumeration algorithm.
#
# The enumeration algorithm is based around growing a seed. A seed
# contains the current subgraph atoms and bonds as well as an exclusion
# set of bonds which cannot be used for future grown. The initial seed
# is the first bond in the fragment, which may potentially grow to use
# the entire fragment. The second seed is the second bond in the
# fragment, which is excluded from using the first bond in future
# growth. The third seed starts from the third bond, which may not use
# the first or second bonds during growth, and so on.
#
#
# A seed can grow along bonds connected to an atom in the seed but which
# aren't already in the seed and aren't in the set of excluded bonds for
# the seed. If there are no such bonds then subgraph enumeration ends
# for this fragment. Given N bonds there are 2**N-1 possible ways to
# grow, which is just the powerset of the available bonds, excluding the
# no-growth case.
#
# This breadth-first growth takes into account all possibilities of using
# the available N bonds so all of those bonds are added to the exclusion
# set of the newly expanded subgraphs.
#
# For performance reasons, the bonds used for growth are separated into
# 'internal' bonds, which connect two atoms already in the subgraph, and
# 'external' bonds, which lead outwards to an atom not already in the
# subgraph.
#
# Each seed growth can add from 0 to N new atoms and bonds. The goal is
# to maximize the subgraph size so the seeds are stored in a priority
# queue, ranked so the seed with the most bonds is processed first. This
# turns the enumeration into something more like a depth-first search.
#
#
# == Prune seeds which aren't found in all of the structures ==
#
# At each stage of seed growth I check that the new seed exists in all
# of the original structures. (Well, all except the one which I
# enumerate over in the first place; by definition that one will match.)
# If it doesn't match then there's no reason to include this seed or any
# larger seeds made from it.
#
# The check is easy; I turn the subgraph into its corresponding SMARTS
# string and use RDKit's normal SMARTS matcher to test for a match.
#
# There are three ways to generate a SMARTS string: 1) arbitrary, 2)
# canonical, 3) hybrid.
#
# I have not tested #1. During most of the development I assumed that
# SMARTS matches across a few hundred structures would be slow, so that
# the best solution is to generate a *canonical* SMARTS and cache the
# match information.
#
# Well, it turns out that my canonical SMARTS match code takes up most
# of the MCS run-time. If I drop the canonicalization step then the
# code averages about 5-10% faster. This isn't the same as #1 - I still
# do the initial atom assignment based on its neighborhood, which is
# like a circular fingerprint of size 2 and *usually* gives a consistent
# SMARTS pattern, which I can then cache.
#
# However, there are times when the non-canonical SMARTS code is slower.
# Obviously one is if there are a lot of structures, and another if is
# there is a lot of symmetry. I'm still working on characterizing this.
#
#
# == Maximize atoms? or bonds? ==
#
# The above algorithm enumerates all subgraphs of the query and
# identifies those subgraphs which are common to all input structures.
#
# It's trivial then to keep track of the current "best" subgraph, which
# can defined as having the subgraph with the most atoms, or the most
# bonds. Both of those options are implemented.
#
# It would not be hard to keep track of all other subgraphs which are
# the same size.
#
# == complete_ring_only implementation ==
#
# The "complete ring only" option is implemented by first enabling the
# "ring-matches-ring-only" option, as otherwise it doesn't make sense.
#
# Second, in order to be a "best" subgraph, all bonds in the subgraph
# which are ring bonds in the original molecule must also be in a ring
# in the subgraph. This is handled as a post-processing step.
#
# (Note: some possible optimizations, like removing ring bonds from
# structure fragments which are not in a ring, are not yet implemented.)
#
#
# == Prune seeds which have no potential for growing large enough  ==
#
# Given a seed, its set of edges available for growth, and the set of
# excluded bonds, figure out the maximum possible growth for the seed.
# If this maximum possible is less than the current best subgraph then
# prune.
#
# This requires a graph search, currently done in Python, which is a bit
# expensive. To speed things up, I precompute some edge information.
# That is, if I know that a given bond is a chain bond (not in a ring)
# then I can calculate the maximum number of atoms and bonds for seed
# growth along that bond, in either direction. However, precomputation
# doesn't take into account the excluded bonds, so after a while the
# predicted value is too high.
#
# Again, I'm still working on characterizing this, and an implementation
# in C++ would have different tradeoffs.

__all__ = ["FindMCS"]

########## Main driver for the MCS code


class MCSResult(object):

  def __init__(self, obj):
    self.numAtoms = obj.num_atoms
    self.numBonds = obj.num_bonds
    self.smarts = obj.smarts
    self.completed = obj.completed

  def __nonzero__(self):
    return self.smarts is not None

  def __repr__(self):
    return "MCSResult(numAtoms=%d, numBonds=%d, smarts=%r, completed=%d)" % (
      self.numAtoms, self.numBonds, self.smarts, self.completed)

  def __str__(self):
    msg = "MCS %r has %d atoms and %d bonds" % (self.smarts, self.numAtoms, self.numBonds)
    if not self.completed:
      msg += " (timed out)"
    return msg


def FindMCS(
  mols,
  minNumAtoms=2,
  maximize=Default.maximize,
  atomCompare=Default.atomCompare,
  bondCompare=Default.bondCompare,
  matchValences=Default.matchValences,
  ringMatchesRingOnly=False,
  completeRingsOnly=False,
  timeout=Default.timeout,
  threshold=None,
):
  """Find the maximum common substructure of a set of molecules

    ***************************************************
    NB: rdkit.Chem.MCS module is deprecated; please use
    rdkit.Chem.rdFMCS instead.
    ***************************************************

    In the simplest case, pass in a list of molecules and get back
    an MCSResult object which describes the MCS:

    >>> from rdkit import Chem
    >>> mols = [Chem.MolFromSmiles("C#CCP"), Chem.MolFromSmiles("C=CCO")]
    >>> from rdkit.Chem import MCS
    >>> MCS.FindMCS(mols)
    MCSResult(numAtoms=2, numBonds=1, smarts='[#6]-[#6]', completed=1)

    The SMARTS '[#6]-[#6]' matches the largest common substructure of
    the input structures. It has 2 atoms and 1 bond. If there is no
    MCS which is at least `minNumAtoms` in size then the result will set
    numAtoms and numBonds to -1 and set smarts to None.

    By default, two atoms match if they are the same element and two
    bonds match if they have the same bond type. Specify `atomCompare`
    and `bondCompare` to use different comparison functions, as in:

    >>> MCS.FindMCS(mols, atomCompare="any")
    MCSResult(numAtoms=3, numBonds=2, smarts='[*]-[*]-[*]', completed=1)
    >>> MCS.FindMCS(mols, bondCompare="any")
    MCSResult(numAtoms=3, numBonds=2, smarts='[#6]~[#6]~[#6]', completed=1)

    An atomCompare of "any" says that any atom matches any other atom,
    "elements" compares by element type, and "isotopes" matches based on
    the isotope label. Isotope labels can be used to implement user-defined
    atom types. A bondCompare of "any" says that any bond matches any
    other bond, and "bondtypes" says bonds are equivalent if and only if
    they have the same bond type.

    A substructure has both atoms and bonds. The default `maximize`
    setting of "atoms" finds a common substructure with the most number
    of atoms. Use maximize="bonds" to maximize the number of bonds.
    Maximizing the number of bonds tends to maximize the number of rings,
    although two small rings may have fewer bonds than one large ring.

    You might not want a 3-valent nitrogen to match one which is 5-valent.
    The default `matchValences` value of False ignores valence information.
    When True, the atomCompare setting is modified to also require that
    the two atoms have the same valency.

    >>> MCS.FindMCS(mols, matchValences=True)
    MCSResult(numAtoms=2, numBonds=1, smarts='[#6v4]-[#6v4]', completed=1)

    It can be strange to see a linear carbon chain match a carbon ring,
    which is what the `ringMatchesRingOnly` default of False does. If
    you set it to True then ring bonds will only match ring bonds.

    >>> mols = [Chem.MolFromSmiles("C1CCC1CCC"), Chem.MolFromSmiles("C1CCCCCC1")]
    >>> MCS.FindMCS(mols)
    MCSResult(numAtoms=7, numBonds=6, smarts='[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]', completed=1)
    >>> MCS.FindMCS(mols, ringMatchesRingOnly=True)
    MCSResult(numAtoms=4, numBonds=3, smarts='[#6](-@[#6])-@[#6]-@[#6]', completed=1)

    You can further restrict things and require that partial rings
    (as in this case) are not allowed. That is, if an atom is part of
    the MCS and the atom is in a ring of the entire molecule then
    that atom is also in a ring of the MCS. Set `completeRingsOnly`
    to True to toggle this requirement and also sets ringMatchesRingOnly
    to True.

    >>> mols = [Chem.MolFromSmiles("CCC1CC2C1CN2"), Chem.MolFromSmiles("C1CC2C1CC2")]
    >>> MCS.FindMCS(mols)
    MCSResult(numAtoms=6, numBonds=6, smarts='[#6]-1-[#6]-[#6](-[#6])-[#6]-1-[#6]', completed=1)
    >>> MCS.FindMCS(mols, ringMatchesRingOnly=True)
    MCSResult(numAtoms=5, numBonds=5, smarts='[#6]-@1-@[#6]-@[#6](-@[#6])-@[#6]-@1', completed=1)
    >>> MCS.FindMCS(mols, completeRingsOnly=True)
    MCSResult(numAtoms=4, numBonds=4, smarts='[#6]-@1-@[#6]-@[#6]-@[#6]-@1', completed=1)

    The MCS algorithm will exhaustively search for a maximum common substructure.
    Typically this takes a fraction of a second, but for some comparisons this
    can take minutes or longer. Use the `timeout` parameter to stop the search
    after the given number of seconds (wall-clock seconds, not CPU seconds) and
    return the best match found in that time. If timeout is reached then the
    `completed` property of the MCSResult will be 0 instead of 1.

    >>> mols = [Chem.MolFromSmiles("Nc1ccccc1"*100), Chem.MolFromSmiles("Nc1ccccccccc1"*100)]
    >>> MCS.FindMCS(mols, timeout=0.1)
    MCSResult(..., completed=0)

    (The MCS after 50 seconds contained 511 atoms.)
    """
  warnings.warn("The rdkit.Chem.MCS module is deprecated; please use rdkit.Chem.rdFMCS instead.",
                DeprecationWarning, stacklevel=2)

  ores = fmcs.fmcs(
    mols,
    minNumAtoms=minNumAtoms,
    maximize=maximize,
    atomCompare=atomCompare,
    bondCompare=bondCompare,
    threshold=threshold,
    matchValences=matchValences,
    ringMatchesRingOnly=ringMatchesRingOnly,
    completeRingsOnly=completeRingsOnly,
    timeout=timeout,
  )
  return MCSResult(ores)


#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest
  import sys
  return doctest.testmod(sys.modules["__main__"],
                         optionflags=doctest.ELLIPSIS + doctest.NORMALIZE_WHITESPACE)


if __name__ == '__main__':
  import sys
  failed, tried = _test()
  sys.exit(failed)
