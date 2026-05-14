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
results (and execution time) it's important to keep an eye on this.)DOC")
      .def_rw("randomSeed", &StereoEnumerationOptions::randomSeed,
              "Seed for random number generator.  Default=-1 means no seed.");

  nb::class_<StereoisomerEnumerator>(m, "StereoisomerEnumerator",
                                     "Stereoisomer enumerator.")
      .def("__init__",
           [](StereoisomerEnumerator *self, const ROMol &mol, bool verbose) {
             new (self)
                 StereoisomerEnumerator(mol, StereoEnumerationOptions(), verbose);
           },
           "mol"_a, "verbose"_a = false)
      .def("__init__",
           [](StereoisomerEnumerator *self, const ROMol &mol,
              const StereoEnumerationOptions &options, bool verbose) {
             new (self) StereoisomerEnumerator(mol, options, verbose);
           },
           "mol"_a, "options"_a, "verbose"_a = false)
      .def("next", &StereoisomerEnumerator::next,
           "Get next isomer in the sequence, or None if at the end.")
      .def("GetStereoisomerCount", &StereoisomerEnumerator::getStereoisomerCount,
           "Get the number of stereoisomers.");
}
