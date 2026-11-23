//
// Copyright (C) 2025-2026 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <nanobind/stl/unique_ptr.h>

#include <GraphMol/EnumerateStereoisomers/EnumerateStereoisomers.h>

namespace nb = nanobind;
using namespace nb::literals;
using namespace RDKit;
using namespace RDKit::EnumerateStereoisomers;

NB_MODULE(rdEnumerateStereoisomers, m) {
  m.doc() =
      R"DOC(Module containing functions to enumerate stereoisomers of a molecule.
Chiral centers and double bonds will be enumerated if unassigned, or,
if the appropriate option is set, if assigned.  Atropisomers will only
be enumerated if assigned.  There is, as yet, no means of finding
unassigned atropisomers.)DOC";

  nb::class_<StereoEnumerationOptions>(m, "StereoEnumerationOptions",
                                       "EnumerateSteroisomers options.")
      .def(nb::init<>())
      .def(nb::init<bool, bool, bool, bool, std::uint64_t, int, bool>(),
           "tryEmbedding"_a = false, "onlyUnassigned"_a = true,
           "onlyStereoGroups"_a = false, "unique"_a = true,
           "maxIsomers"_a = 16384, "randomSeed"_a = 0xC0FFEE,
           "useRingSystemFilter"_a = true)
      .def_rw(
          "tryEmbedding", &StereoEnumerationOptions::tryEmbedding,
          R"DOC(If true, the process attempts to generate a standard RDKit distance geometry
conformation for the stereoisomer.  If this fails, we assume that the stereoisomer is
non-physical and don't return it.  NOTE that this is computationally expensive and is
just a heuristic that could result in stereoisomers being lost.  Default=False)DOC")
      .def_rw(
          "onlyUnassigned", &StereoEnumerationOptions::onlyUnassigned,
          R"DOC(If true, stereocenters which have a specified stereochemistry will not be
perturbed unless they are part of a relative stereo group.  Default=True.)DOC")
      .def_rw(
          "onlyStereoGroups", &StereoEnumerationOptions::onlyStereoGroups,
          R"DOC(If true, only find stereoisomers that differ at the StereoGroups associated with
the molecule.  Default=False.)DOC")
      .def_rw(
          "unique", &StereoEnumerationOptions::unique,
          R"DOC(If true, only stereoisomers that differ in canonical CXSmiles will be
returned.  Default=True.)DOC")
      .def_rw(
          "maxIsomers", &StereoEnumerationOptions::maxIsomers,
          R"DOC(The maximum number of isomers to yield.  If the number of possible isomers
is greater than maxIsomers, a random subset will be yielded.  If 0, there
is no maximum.  Since every additional stereocenter doubles the number of
results (and execution time) it's important to keep an eye on this.
The default value is 16384.")DOC")
      .def_rw(
          "randomSeed", &StereoEnumerationOptions::randomSeed,
          "Seed for random number generator. -1 means no seed. Default=0xC0FFEE")
      .def_rw(
          "useRingSystemFilter", &StereoEnumerationOptions::useRingSystemFilter,
          "filter out stereoisomers containing ring systems that are known to be"
          " impossible in 3D space.  Default=True.");

  nb::class_<StereoisomerEnumerator>(
      m, "StereoisomerEnumerator",
      R"DOC(returns a generator that yields possible stereoisomers for a molecule

    Arguments:
      - m: the molecule to work with
      - options: parameters controlling the enumeration

    If m has stereogroups, they will be expanded

    A small example with 3 chiral atoms and 1 chiral bond (16 theoretical stereoisomers):

    >>> from rdkit import Chem
    >>> from rdkit.Chem.rdEnumerateStereoisomers import StereoisomerEnumerator, StereoEnumerationOptions
    >>> m = Chem.MolFromSmiles('BrC=CC1OC(C2)(F)C2(Cl)C1')
    >>> isomers = tuple(StereoisomerEnumerator(m))
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
    >>> isomers = tuple(StereoisomerEnumerator(m, options=opts))
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
    >>> isomers = tuple(StereoisomerEnumerator(m, options=opts))
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
    >>> isomers = tuple(StereoisomerEnumerator(m))
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
    >>> isomers = tuple(StereoisomerEnumerator(m, options=opts))
    >>> len(isomers)
    16

    Since the result is a generator, we can allow exploring at least parts of very
    large result sets:

    >>> m = Chem.MolFromSmiles('Br' + '[CH](Cl)' * 20 + 'F')
    >>> opts = StereoEnumerationOptions(maxIsomers=0)
    >>> isomers = StereoisomerEnumerator(m, options=opts)
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
    >>> opts = StereoEnumerationOptions(maxIsomers=3,randomSeed=0xf00d)
    >>> isomers = StereoisomerEnumerator(m, options=opts)
    >>> for smi in isomers: #sorted(Chem.MolToSmiles(x, isomericSmiles=True) for x in isomers):
    ...     print(Chem.MolToSmiles(smi))
    F[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)Br
    F[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)Br
    F[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)Br
  )DOC")
      .def(
          "__init__",
          [](StereoisomerEnumerator *self, const ROMol &mol, bool verbose) {
            new (self) StereoisomerEnumerator(mol, StereoEnumerationOptions(),
                                              verbose);
          },
          "mol"_a, "verbose"_a = false)
      .def(
          "__init__",
          [](StereoisomerEnumerator *self, const ROMol &mol,
             const StereoEnumerationOptions &options, bool verbose) {
            new (self) StereoisomerEnumerator(mol, options, verbose);
          },
          "mol"_a, "options"_a, "verbose"_a = false)
      .def("__iter__", [](nb::object self) { return self; })
      .def(
          "__next__",
          [](StereoisomerEnumerator &self) {
            auto isomer = self.next();
            if (!isomer) {
              throw nb::stop_iteration();
            }
            return isomer;
          },
          "Get the next isomer in the sequence. Raises StopIteration at the end.")
      .def("next", &StereoisomerEnumerator::next,
           "Get next isomer in the sequence, or None if at the end.")
      .def("GetStereoisomerCount",
           &StereoisomerEnumerator::getStereoisomerCount,
           "Get the number of stereoisomers.");
}
