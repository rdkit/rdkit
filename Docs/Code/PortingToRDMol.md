# Porting code to use RDMol

`RDMol` is the new data structure for representing molecule data in the RDKit.
However, to maintain as much compatibility as possible, `RWMol`, `ROMol`,
`Atom`, and `Bond` are now wrappers that still provide compatibility interfaces
that closely mimic the behaviour of the previous data structures, while their
data is now primarily managed by `RDMol`.

The new data structure is more memory efficient and incurs less memory
management overhead in many situations. The compatibility interface should
still function as it did before in a vast majority of situations, however,
the compatibility interface incurs a significant performance penalty in some
cases, and the new interface can allow for some optimizations that weren't
available before, so it's recommended to port existing code to use `RDMol`.

This guide provides common examples of how to port existing code. Code can
be ported in smaller chunks by using `RDMol::asRWMol()` or `asROMol()` to get
an `RWMol` or `ROMol` wrapper from an `RDMol`, and using `ROMol::asRDMol()`
to get the `RDMol` from an `ROMol`. This allows porting of performance-critical
or high priority code first, while leaving porting of other code until later.
Care does need to be taken to not access both interfaces at the same time,
whether via multi-threading, or via caching pointers from one interface that
may or may not become stale upon using the other interface. Apart from that,
switching between the interfaces is simple. Supposing a function `foo(ROMol&)`
has been ported to `foo(RDMol&)`, but the calling code hasn't been ported yet,
the calling code can be updated from using `foo(romol)` to using
`foo(romol.asRDMol())`. Similarly, if the calling code has been ported, but the
function hasn't, the call to `foo` would look like `foo(rdmol.asROMol())`.

## Accessing atom and bond data

Basic atom and bond data that were previously accessed via `Atom` and `Bond`
objects are now stored in vectors of `AtomData` and `BondData`, managed by
and accessed via `RDMol`.  For example,

```
const Bond *bond = romol.getBondWithIdx(bondIdx);
uint32_t beginAtomIdx = bond->getBeginAtomIdx();
uint32_t endAtomIdx = bond->getEndAtomIdx();
Bond::BondType type = bond->getBondType();
const Atom *atom = romol.getAtomWithIdx(beginAtomIdx);
int atomicNum = atom->getAtomicNum();
```

becomes

```
const BondData &bond = rdmol.getBond(bondIdx);
uint32_t beginAtomIdx = bond.getBeginAtomIdx();
uint32_t endAtomIdx = bond.getEndAtomIdx();
BondEnums::BondType type = bond.getBondType();
const AtomData &atom = rdmol.getAtom(beginAtomIdx);
uint32_t atomicNum = atom.getAtomicNum();
```

Accessing a bond by index is now constant time, instead of a linear search, so
simple iteration over bonds by index is efficient now. Iteration with range-
based for loops is still supported, but for performance critical code, using
indices directly may be slightly faster. For example:

```
for (Atom *atom : romol.atoms()) {
  int atomicNum = atom->getAtomicNum();
  ...
}
for (Bond *bond : romol.bonds()) {
  uint32_t beginAtomIdx = bond->getBeginAtomIdx();
  uint32_t endAtomIdx = bond->getEndAtomIdx();
  ...
}
```

can become

```
for (RDMolAtom atom : rdmol.atoms()) {
  int atomicNum = atom.data().getAtomicNum();
  ...
}
for (RDMolBond bond : rdmol.bonds()) {
  uint32_t beginAtomIdx = bond.data().getBeginAtomIdx();
  uint32_t endAtomIdx = bond.data().getEndAtomIdx();
  ...
}
```

or with a standard for loop

```
for (uint32_t atomIdx = 0, numAtoms = rdmol.getNumAtoms();
    atomIdx < numAtoms; ++atomIdx) {
  const AtomData &atom = rdmol.getAtom(atomIdx);
  int atomicNum = atom.getAtomicNum();
  ...
}
for (uint32_t bondIdx = 0, numBonds = rdmol.getNumBonds();
    bondIdx < numBonds; ++bondIdx) {
  const BondData &bond = rdmol.getBond(bondIdx);
  uint32_t beginAtomIdx = bond.data().getBeginAtomIdx();
  uint32_t endAtomIdx = bond.data().getEndAtomIdx();
  ...
}
```

The new interface no longer allows getting mutable access to atom and bond
data from a `const RDMol`, so for the `const ` case, the range-based for loops
above would have `ConstRDMolAtom` or `ConstRDMolBond` in place of `RDMolAtom`
or `RDMolBond`. These 4 classes are simple wrappers of a pointer to `RDMol` and
an index of either an atom or a bond, to be able to passed around fairly lightly
in places where `Atom` or `Bond` pointers were used before. They are invalid
if default constructed, which acts like a null pointer, and this can be checked
by calling `isValid()` on them. A reference to the `RDMol` is returned by
`mol()`, and the atom or bond index is returned by `index()`. In situations
where only one molecule is involved, using only an index is usually more
efficient than using these wrappers. Unlike `Atom` and `Bond` pointers, indices
and these wrappers are not persistent after removing an atom or bond, so be
careful not to use out of date indices in those cases. Similarly, don't cache
pointers to `AtomData` or `BondData` if atoms or bonds might be added or
removed, because the vector containing that data may be reallocated,
invalidating those pointers.


Bond stereo atoms used to be a `std::vector<int>` on each `Bond`, but only 0 or
2 atoms were valid sizes, so this is now represented as a fixed array of 2 atom
indices, with -1 indicating that the indices are unused. However, in order to
support the previous format in the compatibility interface, the new
representation is accessed via `RDMol::setBondStereoAtoms`, instead of accessing
the indices in `BondData` directly, in order to synchronize between the
compatibility data and the new data, as needed. For example,
`bond->setStereoAtoms(beginIdx, endIdx);` can become
`rdmol.setBondStereoAtoms(bondIdx, beginIdx, endIdx);` Because only sizes 0 or 2
are valid, code that previously used `bond->getStereoAtoms().clear();` can now
call `rdmol.clearBondStereoAtoms();` and code that used `getStereoAtoms()` to
set 2 atom indices can be adapted to call `RDMol::setBondStereoAtoms`.
Similarly, there are now `RDMol::hasBondStereoAtoms(bondIdx)` and
`RDMol::getBondStereoAtoms(bondIdx)` for reading the bond stereo atoms.
`RDMol::getBondStereoAtoms(bondIdx)` currently always returns non-null, but
the atom indices will be invalid (`atomindex_t(-1)`) if the bond does not have
stereo atoms.

## Neighbours and graph structure

Previously, the graph representing the connections between atoms used the Boost
Graph library, which limited access to using iterators provided by that library.
Similar range-based for loops are still possible with the new graph structure,
represented using flattened vectors, but now some further optimizations are
possible for performance-critical code. For example, iterating over the bonds
of an atom:

```
for (Bond *bond : romol.atomBonds(atom)) {
  ...
}
```

can become the following, using the atom index, instead of a pointer,

```
for (RDMolBond bond : rdmol.atomBonds(atomIdx)) {
  ...
}
```

or, for more direct access to the bond indices

```
// The "auto" here represents "const uint32_t*"
auto [begin, end] = rdmol.getAtomBonds(atomIdx);
for (; begin != end; ++begin) {
  uint32_t bondIdx = *begin;
  const BondData &bond = rdmol.getBond(bondIdx);
  ...
}
```

For iterating over the neighbours of an atom, instead of the bonds, the
corresponding functions are `atomNeighbours` or `getAtomNeighbours`.
To get the number of bonds/neighbours of an atom, `rdmol.getAtomDegree(atomIdx)`
replaces `romol.getAtomDegree(atom)` or `atom->getDegree()`. Like the atom
and bond iterators in the previous section, `const RDMol` no longer allows
acquiring mutable access to atom and bond data, so in the `const` case,
`atomBonds` and `atomNeighbours` will yield `ConstRDMolBond` or
`ConstRDMolAtom`.

## Constructing a molecule

Just like `ROMol` or `RWMol`, `RDMol` can be default constructed to make an
empty molecule. However, in the new interface, atom and bond data are managed by
`RDMol`, instead of being managed by separately allocated `Atom` and `Bond`
objects, so for example

```
RWMol rwmol;
for (uint32_t atomIdx = 0; atomIdx < numAtoms; ++atomIdx) {
  Atom *atom = new Atom(atomicNums[atomIdx]);
  atom->setFormalCharge(atomCharges[atomIdx]);
  rwmol.addAtom(atom, true, /*takeOwnership=*/true);
}
for (uint32_t bondIdx = 0; bondIdx < numBonds; ++bondIdx) {
  Bond *bond = new Bond(bondTypes[bondIdx]);
  bond->setBeginAtomIdx(bondAtoms[bondIdx].first);
  bond->setEndAtomIdx(bondAtoms[bondIdx].second);
  bond->setStereo(bondStereos[bondIdx]);
  rwmol.addBond(bond, /*takeOwnership=*/true);
}
```

could become

```
RDMol rdmol;
for (uint32_t atomIdx = 0; atomIdx < numAtoms; ++atomIdx) {
  AtomData &atom = rdmol.addAtom();
  atom.setAtomicNum(atomicNums[atomIdx]);
  atom.setFormalCharge(atomCharges[atomIdx]);
}
for (uint32_t bondIdx = 0; bondIdx < numBonds; ++bondIdx) {
  BondData &bond = rdmol.addBond(
    bondAtoms[bondIdx].first, bondAtoms[bondIdx].second, bondTypes[bondIdx]);
  bond.setStereo(bondStereos[bondIdx]);
}
```

## Properties

In most situations, properties, "props", are not needed. If you can use a
separate structure to store the relevant data outside of the molecule, it
data structure, except for a few cases, such as using a property to keep track
will likely be more efficient than using properties. However, properties are
used in many places, and `RDMol` provides some new features to make the most
common uses of them more efficient.

Because `RDMol` now manages atom and bond data, the common case of a property
being added to all properties can be handled more efficiently, using an array,
especially for simple types `bool`, `char`, `int32_t`, `uint32t`, `int64_t`,
`uint64_t`, `float`, and `double`. If a string property has a limited set of
valid strings, it's recommended to replace the strings with an enum, so that
an integer property can be used instead. Properties can still be present or
absent from each atom or bond, but there's less overhead now for each additional
atom or bond having a property after the first one within a molecule. That
means there's no advantage to having a `bool` property on certain atoms and
not others, over having the property on all atoms.

The other large change with properties is that, because property names are
frequently known at compile time, a new shared string `PropToken` class is now
used instead of using `std::string` for property names, to reduce string
copying. For faster comparison, it also caches a hash of the string.

Looking at an example:

```

```



## Ring info

With the new interface, ring info is now stored in a flattened representation,
`RingInfoCache`, which is especially efficient for the case of multiple rings.
`ringBegins` contains the indices into `atomsInRings` and `bondsInRings` where
each ring begins, so that `atomsInRings` and `bondsInRings` can be single
vectors for containing data for all rings. This changes how to access the ring
data. For example:

```
const RingInfo *info = romol.getRingInfo();
const uint32_t numRings = info->numRings();
const std::vector<std::vector<int>> &atomRings = info->atomRings();
const std::vector<std::vector<int>> &bondRings = info->bondRings();
for (uint32_t ringIdx = 0; ringIdx < numRings; ++ringIdx) {
  const std::vector<int> &atomRing = atomRings[ringIdx];
  const std::vector<int> &bondRing = bondRings[ringIdx];
  const uint32_t ringSize = atomRing.size();
  for (int atomIdx : atomRing) {
    ...
  }
  for (int bondIdx : bondRing) {
    ...
  }
}
```

could be ported as:

```
const RingInfoCache &info = rdmol.getRingInfo();
const uint32_t numRings = info.numRings();
const std::vector<uint32_t> &atomRings = info.atomsInRings;
const std::vector<uint32_t> &bondRings = info.bondsInRings;
for (uint32_t ringIdx = 0; ringIdx < numRings; ++ringIdx) {
  const uint32_t ringBegin = info.ringBegins[ringIdx];
  const uint32_t ringEnd = info.ringBegins[ringIdx + 1];
  const uint32_t ringSize = ringEnd - ringBegin;
  for (uint32_t idx = ringBegin; idx != ringEnd; ++idx) {
    const uint32_t atomIdx = atomRings[idx];
    ...
  }
  for (uint32_t idx = ringBegin; idx != ringEnd; ++idx) {
    const uint32_t bondIdx = bondRings[idx];
    ...
  }
}
```

`ringBegins` has an extra value at the end, to avoid needing a bounds check
for `ringIdx + 1`. For iterating over the rings of a particular atom or the
rings of a particular bond, `atomMembershipBegins` indicates where each atom's
list of rings begins in `atomMemberships`, and `bondMembershipBegins` indicates
where each bond's list of rings begins in `bondMemberships`. For example:

```
for (const Atom *atom : romol.atoms()) {
  const uint32_t atomIdx = atom->getIdx();
  const std::vector<int> &rings = info->atomMembers(atomIdx);
  const atomNumRings = rings.size();
  for (int ringIdx : rings) {
    // Visit all bonds in rings that contain atom
    const std::vector<int> &bondRing = info->bondRings()[ringIdx];
    for (int bondIdx : bondRing) {
      ...
    }
  }
}
```

could be ported to:

```
for (const uint32_t atomIdx = 0, numAtoms = rdmol.getNumAtoms();
    atomIdx < numAtoms; ++atomIdx) {
  const uint32_t ringsBegin = info.atomMembershipBegins[atomIdx];
  const uint32_t ringsEnd = info.atomMembershipBegins[atomIdx + 1];
  const atomNumRings = ringsEnd - ringsBegin;
  for (uint32_t idx = ringsBegin; idx != ringsEnd; ++idx) {
    const uint32_t ringIdx = info.atomMemberships[idx];
    // Visit all bonds in rings that contain atom
    const uint32_t ringBegin = info.ringBegins[ringIdx];
    const uint32_t ringEnd = info.ringBegins[ringIdx + 1];
    for (uint32_t idx2 = ringBegin; idx2 != ringEnd; ++idx2) {
      const uint32_t bondIdx = info.bondRings[idx2];
      ...
    }
  }
}
```

Like `ringBegins`, `atomMembershipBegins` has an extra value at the end, to
avoid needing a bounds check for `atomIdx + 1`.


## Conformer data

## Queries

## Safety of multi-threading

## Known compatibility issues

