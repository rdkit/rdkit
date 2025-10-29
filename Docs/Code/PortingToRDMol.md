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
const Atom *atom = romol.getAtomWithIdx(atomIdx);
int atomicNum = atom->getAtomicNum();
```

becomes

```
const AtomData &atom = rdmol.getAtom(atomIdx);
uint32_t atomicNum = atom.getAtomicNum();
```

## Constructing a molecule

## Neighbours and graph structure

## Properties

## Iteration

## Safety of multi-threading

