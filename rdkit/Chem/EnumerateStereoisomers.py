import random

from rdkit.Chem import rdEnumerateStereoisomers
from rdkit import Chem

MAX_RAND = 2**31-1
class StereoEnumerationOptions(object):
  """
          - tryEmbedding: if set the process attempts to generate a standard RDKit distance geometry
            conformation for the stereisomer. If this fails, we assume that the stereoisomer is
            non-physical and don't return it. NOTE that this is computationally expensive and is
            just a heuristic that could result in stereoisomers being lost.

          - onlyUnassigned: if set (the default), stereocenters which have specified stereochemistry
            will not be perturbed unless they are part of a relative stereo
            group.

          - maxIsomers: the maximum number of isomers to yield, if the
            number of possible isomers is greater than maxIsomers, a
            random subset will be yielded. If 0, all isomers are
            yielded. Since every additional stereo center doubles the
            number of results (and execution time) it's important to
            keep an eye on this.

          - onlyStereoGroups: Only find stereoisomers that differ at the
            StereoGroups associated with the molecule.
    """
  __slots__ = ('opts',)

  def __init__(self, tryEmbedding=False, onlyUnassigned=True, maxIsomers=1024, rand=None,
               unique=True, onlyStereoGroups=False):
    self.opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    self.opts.tryEmbedding = tryEmbedding
    self.opts.onlyUnassigned = onlyUnassigned
    self.opts.onlyStereoGroups = onlyStereoGroups
    self.opts.maxIsomers = maxIsomers
    if rand is not None:
      if int(rand) > MAX_RAND:
        rand = MAX_RAND
      self.opts.randomSeed = rand
    self.opts.unique = unique


def GetStereoisomerCount(m, options=StereoEnumerationOptions()):
  """ returns an estimate (upper bound) of the number of possible stereoisomers for a molecule

   Arguments:
      - m: the molecule to work with
      - options: parameters controlling the enumeration


    >>> from rdkit import Chem
    >>> from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
    >>> m = Chem.MolFromSmiles('BrC(Cl)(F)CCC(O)C')
    >>> GetStereoisomerCount(m)
    4
    >>> m = Chem.MolFromSmiles('CC(Cl)(O)C')
    >>> GetStereoisomerCount(m)
    1

    double bond stereochemistry is also included:

    >>> m = Chem.MolFromSmiles('BrC(Cl)(F)C=CC(O)C')
    >>> GetStereoisomerCount(m)
    8

    """
  enumerator = rdEnumerateStereoisomers.StereoisomerEnumerator(m, options.opts)
  return enumerator.GetStereoisomerCount()


def EnumerateStereoisomers(m, options=StereoEnumerationOptions(), verbose=False):
  r""" returns a generator that yields possible stereoisomers for a molecule

    Arguments:
      - m: the molecule to work with
      - options: parameters controlling the enumeration
      - verbose: toggles how verbose the output is

    If m has stereogroups, they will be expanded

    A small example with 3 chiral atoms and 1 chiral bond (16 theoretical stereoisomers):

    >>> from rdkit import Chem
    >>> from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
    >>> m = Chem.MolFromSmiles('BrC=CC1OC(C2)(F)C2(Cl)C1')
    >>> isomers = tuple(EnumerateStereoisomers(m))
    >>> len(isomers)
    16
    >>> for smi in sorted(Chem.MolToSmiles(x, isomericSmiles=True) for x in isomers):
    ...     print(smi)
    ...
    F[C@@]12C[C@@]1(Cl)C[C@@H](/C=C/Br)O2
    F[C@@]12C[C@@]1(Cl)C[C@@H](/C=C\Br)O2
    F[C@@]12C[C@@]1(Cl)C[C@H](/C=C/Br)O2
    F[C@@]12C[C@@]1(Cl)C[C@H](/C=C\Br)O2
    F[C@@]12C[C@]1(Cl)C[C@@H](/C=C/Br)O2
    F[C@@]12C[C@]1(Cl)C[C@@H](/C=C\Br)O2
    F[C@@]12C[C@]1(Cl)C[C@H](/C=C/Br)O2
    F[C@@]12C[C@]1(Cl)C[C@H](/C=C\Br)O2
    F[C@]12C[C@@]1(Cl)C[C@@H](/C=C/Br)O2
    F[C@]12C[C@@]1(Cl)C[C@@H](/C=C\Br)O2
    F[C@]12C[C@@]1(Cl)C[C@H](/C=C/Br)O2
    F[C@]12C[C@@]1(Cl)C[C@H](/C=C\Br)O2
    F[C@]12C[C@]1(Cl)C[C@@H](/C=C/Br)O2
    F[C@]12C[C@]1(Cl)C[C@@H](/C=C\Br)O2
    F[C@]12C[C@]1(Cl)C[C@H](/C=C/Br)O2
    F[C@]12C[C@]1(Cl)C[C@H](/C=C\Br)O2

    Because the molecule is constrained, not all of those isomers can
    actually exist. We can check that:

    >>> opts = StereoEnumerationOptions(tryEmbedding=True)
    >>> isomers = tuple(EnumerateStereoisomers(m, options=opts))
    >>> len(isomers)
    8
    >>> for smi in sorted(Chem.MolToSmiles(x,isomericSmiles=True) for x in isomers):
    ...     print(smi)
    ...
    F[C@@]12C[C@]1(Cl)C[C@@H](/C=C/Br)O2
    F[C@@]12C[C@]1(Cl)C[C@@H](/C=C\Br)O2
    F[C@@]12C[C@]1(Cl)C[C@H](/C=C/Br)O2
    F[C@@]12C[C@]1(Cl)C[C@H](/C=C\Br)O2
    F[C@]12C[C@@]1(Cl)C[C@@H](/C=C/Br)O2
    F[C@]12C[C@@]1(Cl)C[C@@H](/C=C\Br)O2
    F[C@]12C[C@@]1(Cl)C[C@H](/C=C/Br)O2
    F[C@]12C[C@@]1(Cl)C[C@H](/C=C\Br)O2

    Or we can force the output to only give us unique isomers:

    >>> m = Chem.MolFromSmiles('FC(Cl)C=CC=CC(F)Cl')
    >>> opts = StereoEnumerationOptions(unique=True)
    >>> isomers = tuple(EnumerateStereoisomers(m, options=opts))
    >>> len(isomers)
    10
    >>> for smi in sorted(Chem.MolToSmiles(x,isomericSmiles=True) for x in isomers):
    ...     print(smi)
    ...
    F[C@@H](Cl)/C=C/C=C/[C@@H](F)Cl
    F[C@@H](Cl)/C=C\C=C/[C@@H](F)Cl
    F[C@@H](Cl)/C=C\C=C\[C@@H](F)Cl
    F[C@H](Cl)/C=C/C=C/[C@@H](F)Cl
    F[C@H](Cl)/C=C/C=C/[C@H](F)Cl
    F[C@H](Cl)/C=C/C=C\[C@@H](F)Cl
    F[C@H](Cl)/C=C\C=C/[C@@H](F)Cl
    F[C@H](Cl)/C=C\C=C/[C@H](F)Cl
    F[C@H](Cl)/C=C\C=C\[C@@H](F)Cl
    F[C@H](Cl)/C=C\C=C\[C@H](F)Cl

    By default the code only expands unspecified stereocenters:

    >>> m = Chem.MolFromSmiles('BrC=C[C@H]1OC(C2)(F)C2(Cl)C1')
    >>> isomers = tuple(EnumerateStereoisomers(m))
    >>> len(isomers)
    8
    >>> for smi in sorted(Chem.MolToSmiles(x,isomericSmiles=True) for x in isomers):
    ...     print(smi)
    ...
    F[C@@]12C[C@@]1(Cl)C[C@@H](/C=C/Br)O2
    F[C@@]12C[C@@]1(Cl)C[C@@H](/C=C\Br)O2
    F[C@@]12C[C@]1(Cl)C[C@@H](/C=C/Br)O2
    F[C@@]12C[C@]1(Cl)C[C@@H](/C=C\Br)O2
    F[C@]12C[C@@]1(Cl)C[C@@H](/C=C/Br)O2
    F[C@]12C[C@@]1(Cl)C[C@@H](/C=C\Br)O2
    F[C@]12C[C@]1(Cl)C[C@@H](/C=C/Br)O2
    F[C@]12C[C@]1(Cl)C[C@@H](/C=C\Br)O2

    But we can change that behavior:

    >>> opts = StereoEnumerationOptions(onlyUnassigned=False)
    >>> isomers = tuple(EnumerateStereoisomers(m, options=opts))
    >>> len(isomers)
    16

    Since the result is a generator, we can allow exploring at least parts of very
    large result sets:

    >>> m = Chem.MolFromSmiles('Br' + '[CH](Cl)' * 20 + 'F')
    >>> opts = StereoEnumerationOptions(maxIsomers=0)
    >>> isomers = EnumerateStereoisomers(m, options=opts)
    >>> for x in range(5):
    ...   print(Chem.MolToSmiles(next(isomers),isomericSmiles=True))  # doctest: +SKIP
    F[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)Br
    F[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)Br
    F[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)Br
    F[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)Br
    F[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)Br

    Or randomly sample a small subset. Note that if we want that sampling to be consistent
    across runs we need to provide a random number seed:

    >>> m = Chem.MolFromSmiles('Br' + '[CH](Cl)' * 20 + 'F')
    >>> opts = StereoEnumerationOptions(maxIsomers=3,rand=0xf00d)
    >>> isomers = EnumerateStereoisomers(m, options=opts)
    >>> for smi in isomers: #sorted(Chem.MolToSmiles(x, isomericSmiles=True) for x in isomers):
    ...     print(Chem.MolToSmiles(smi))
    F[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)Br
    F[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)Br
    F[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)Br

    Additionally, isomers are marked with their isomerSet.  When using enhanced stereo labels,
    STEREO_ORs indicate that there is the possibiity of two batches, one with R and one with S.
    Each isomerSet indicates the isomers that are seperated into batches.  For example, the
    following molecule has one STEREO_OR group so should make two isomerSets:

    >>> m = Chem.MolFromSmiles('C[C@@H](Br)[C@@H](C)[C@H](C)F |wD:3.3,5.5,wU:1.0,a:1,o1:3,5|')
    >>> isomers = [(Chem.MolToSmiles(ism), ism) for ism in EnumerateStereoisomers(m)]
    >>> for smiles, ism in sorted(isomers):
    ...    print("isomerSet:", ism.GetIntProp("isomerSet"), smiles)
    isomerSet: 1 C[C@@H]([C@H](C)F)[C@@H](C)Br
    isomerSet: 0 C[C@H]([C@@H](C)F)[C@@H](C)Br

    """
  enumerator = rdEnumerateStereoisomers.StereoisomerEnumerator(m, options.opts)
  return iter(enumerator)

#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest
  import sys
  return doctest.testmod(sys.modules["__main__"])


if __name__ == '__main__':
  import sys
  failed, tried = _test()
  sys.exit(failed)
