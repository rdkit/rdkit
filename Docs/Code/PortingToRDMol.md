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

## Properties

## Iteration

## Ring info

## Conformer data

## Safety of multi-threading

## Known compatibility issues

