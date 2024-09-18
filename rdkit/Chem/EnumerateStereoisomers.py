import random

from rdkit import Chem
from rdkit.Chem.rdDistGeom import EmbedMolecule


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
  __slots__ = ('tryEmbedding', 'onlyUnassigned', 'onlyStereoGroups', 'maxIsomers', 'rand', 'unique')

  def __init__(self, tryEmbedding=False, onlyUnassigned=True, maxIsomers=1024, rand=None,
               unique=True, onlyStereoGroups=False):
    self.tryEmbedding = tryEmbedding
    self.onlyUnassigned = onlyUnassigned
    self.onlyStereoGroups = onlyStereoGroups
    self.maxIsomers = maxIsomers
    self.rand = rand
    self.unique = unique


class _BondFlipper(object):

  def __init__(self, bond):
    self.bond = bond

  def flip(self, flag):
    if flag:
      self.bond.SetStereo(Chem.BondStereo.STEREOCIS)
    else:
      self.bond.SetStereo(Chem.BondStereo.STEREOTRANS)


class _AtomFlipper(object):

  def __init__(self, atom):
    self.atom = atom

  def flip(self, flag):
    if flag:
      self.atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    else:
      self.atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)


class _StereoGroupFlipper(object):

  def __init__(self, group):
    self._original_parities = [(a, a.GetChiralTag()) for a in group.GetAtoms()]

  def flip(self, flag):
    if flag:
      for a, original_parity in self._original_parities:
        a.SetChiralTag(original_parity)
    else:
      for a, original_parity in self._original_parities:
        if original_parity == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
          a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif original_parity == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
          a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)


def _getFlippers(mol, options):
  sinfo = Chem.FindPotentialStereo(mol)
  flippers = []
  if not options.onlyStereoGroups:
    for si in sinfo:
      if options.onlyUnassigned and si.specified not in (Chem.StereoSpecified.Unspecified,
                                                         Chem.StereoSpecified.Unknown):
        continue
      if si.type == Chem.StereoType.Atom_Tetrahedral:
        flippers.append(_AtomFlipper(mol.GetAtomWithIdx(si.centeredOn)))
      elif si.type == Chem.StereoType.Bond_Double:
        bnd = mol.GetBondWithIdx(si.centeredOn)
        if not bnd.GetStereoAtoms():
          if si.controllingAtoms[0] == Chem.StereoInfo.NOATOM or \
            si.controllingAtoms[2] == Chem.StereoInfo.NOATOM:
            continue
          bnd.SetStereoAtoms(si.controllingAtoms[0], si.controllingAtoms[2])
        flippers.append(_BondFlipper(mol.GetBondWithIdx(si.centeredOn)))
      ## FIX: support atropisomers

  if options.onlyUnassigned:
    # otherwise these will be counted twice
    for group in mol.GetStereoGroups():
      if group.GetGroupType() != Chem.StereoGroupType.STEREO_ABSOLUTE:
        flippers.append(_StereoGroupFlipper(group))

  return flippers


class _RangeBitsGenerator(object):

  def __init__(self, nCenters):
    self.nCenters = nCenters

  def __iter__(self):
    for val in range(2**self.nCenters):
      yield val


class _UniqueRandomBitsGenerator(object):

  def __init__(self, nCenters, maxIsomers, rand):
    self.nCenters = nCenters
    self.maxIsomers = maxIsomers
    self.rand = rand
    self.already_seen = set()

  def __iter__(self):
    # note: important that this is not 'while True' otherwise it
    # would be possible to have an infinite loop caused by all
    # isomers failing the embedding process
    while len(self.already_seen) < 2**self.nCenters:
      bits = self.rand.getrandbits(self.nCenters)
      if bits in self.already_seen:
        continue

      self.already_seen.add(bits)
      yield bits


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
  tm = Chem.Mol(m)
  flippers = _getFlippers(tm, options)
  return 2**len(flippers)


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
    ...   print(Chem.MolToSmiles(next(isomers),isomericSmiles=True))
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
    F[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)Br
    F[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)Br
    F[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)Br

    """
  tm = Chem.Mol(m)
  for atom in tm.GetAtoms():
    atom.ClearProp("_CIPCode")
  for bond in tm.GetBonds():
    if bond.GetBondDir() == Chem.BondDir.EITHERDOUBLE:
      bond.SetBondDir(Chem.BondDir.NONE)
  flippers = _getFlippers(tm, options)
  nCenters = len(flippers)
  if not nCenters:
    yield tm
    return

  if (options.maxIsomers == 0 or 2**nCenters <= options.maxIsomers):
    bitsource = _RangeBitsGenerator(nCenters)
  else:
    if options.rand is None:
      # deterministic random seed invariant to input atom order
      seed = hash(tuple(sorted([(a.GetDegree(), a.GetAtomicNum()) for a in tm.GetAtoms()])))
      rand = random.Random(seed)
    elif isinstance(options.rand, random.Random):
      # other implementations of Python random number generators
      # can inherit from this class to pick up utility methods
      rand = options.rand
    else:
      rand = random.Random(options.rand)

    bitsource = _UniqueRandomBitsGenerator(nCenters, options.maxIsomers, rand)

  isomersSeen = set()
  numIsomers = 0
  for bitflag in bitsource:
    for i in range(nCenters):
      flag = bool(bitflag & (1 << i))
      flippers[i].flip(flag)

    # from this point on we no longer need the stereogroups (if any are there), so
    # remove them:
    if tm.GetStereoGroups():
      isomer = Chem.RWMol(tm)
      isomer.SetStereoGroups([])
    else:
      isomer = Chem.Mol(tm)
    Chem.SetDoubleBondNeighborDirections(isomer)
    isomer.ClearComputedProps(includeRings=False)

    Chem.AssignStereochemistry(isomer, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    if options.unique:
      cansmi = Chem.MolToSmiles(isomer, isomericSmiles=True)
      if cansmi in isomersSeen:
        continue

      isomersSeen.add(cansmi)

    if options.tryEmbedding:
      ntm = Chem.AddHs(isomer)
      # mask bitflag to fit within C++ int.
      cid = EmbedMolecule(ntm, randomSeed=(bitflag & 0x7fffffff))
      if cid >= 0:
        conf = Chem.Conformer(isomer.GetNumAtoms())
        for aid in range(isomer.GetNumAtoms()):
          conf.SetAtomPosition(aid, ntm.GetConformer().GetAtomPosition(aid))
        isomer.AddConformer(conf)
    else:
      cid = 1
    if cid >= 0:
      yield isomer
      numIsomers += 1
      if options.maxIsomers != 0 and numIsomers >= options.maxIsomers:
        break
    elif verbose:
      print("%s    failed to embed" % (Chem.MolToSmiles(isomer, isomericSmiles=True)))


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
