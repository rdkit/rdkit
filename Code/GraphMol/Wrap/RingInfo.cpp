//
//  Copyright (C) Greg Landrum 2007-2017
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define NO_IMPORT_ARRAY
#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>

#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>

namespace python = boost::python;

namespace {
using namespace RDKit;
python::object atomRings(const RingInfo *self) {
  python::list res;
  for (const auto &ring : self->atomRings()) {
    res.append(python::tuple(ring));
  }
  return python::tuple(res);
}
python::object bondRings(const RingInfo *self) {
  python::list res;
  for (const auto &ring : self->bondRings()) {
    res.append(python::tuple(ring));
  }
  return python::tuple(res);
}
python::object atomMembers(const RingInfo *self, unsigned int idx) {
  return python::tuple(self->atomMembers(idx));
}
python::object bondMembers(const RingInfo *self, unsigned int idx) {
  return python::tuple(self->bondMembers(idx));
}
python::object atomRingSizes(const RingInfo *self, unsigned int idx) {
  return python::tuple(self->atomRingSizes(idx));
}
python::object bondRingSizes(const RingInfo *self, unsigned int idx) {
  return python::tuple(self->bondRingSizes(idx));
}

python::object atomRingFamilies(const RingInfo *self) {
  python::list res;
  for (const auto &ring : self->atomRingFamilies()) {
    res.append(python::tuple(ring));
  }
  return python::tuple(res);
}
python::object bondRingFamilies(const RingInfo *self) {
  python::list res;
  for (const auto &ring : self->bondRingFamilies()) {
    res.append(python::tuple(ring));
  }
  return python::tuple(res);
}

void addRing(RingInfo *self, python::object atomRing, python::object bondRing) {
  unsigned int nAts = boost::python::len(atomRing);
  unsigned int nBnds = boost::python::len(bondRing);
  if (nAts != nBnds) {
    throw_value_error("list sizes must match");
  }
  if (!self->isInitialized()) {
    self->initialize();
  }
  INT_VECT aring(nAts);
  INT_VECT bring(nAts);
  for (unsigned int i = 0; i < nAts; ++i) {
    aring[i] = python::extract<int>(atomRing[i])();
    bring[i] = python::extract<int>(bondRing[i])();
  }
  self->addRing(aring, bring);
}
}  // namespace

namespace RDKit {
std::string classDoc = "contains information about a molecule's rings\n";

struct ringinfo_wrapper {
  static void wrap() {
    python::class_<RingInfo>("RingInfo", classDoc.c_str(), python::no_init)
        .def("IsAtomInRingOfSize", &RingInfo::isAtomInRingOfSize,
             python::args("self", "idx", "size"))
        .def("MinAtomRingSize", &RingInfo::minAtomRingSize,
             python::args("self", "idx"))
        .def("AreAtomsInSameRing", &RingInfo::areAtomsInSameRing,
             python::args("self", "idx1", "idx2"))
        .def("AreAtomsInSameRingOfSize", &RingInfo::areAtomsInSameRingOfSize,
             python::args("self", "idx1", "idx2", "size"))
        .def("IsBondInRingOfSize", &RingInfo::isBondInRingOfSize,
             python::args("self", "idx", "size"))
        .def("MinBondRingSize", &RingInfo::minBondRingSize,
             python::args("self", "idx"))
        .def("AreBondsInSameRing", &RingInfo::areBondsInSameRing,
             python::args("self", "idx1", "idx2"))
        .def("AreBondsInSameRingOfSize", &RingInfo::areBondsInSameRingOfSize,
             python::args("self", "idx1", "idx2", "size"))
        .def("NumAtomRings", &RingInfo::numAtomRings,
             python::args("self", "idx"))
        .def("NumBondRings", &RingInfo::numBondRings,
             python::args("self", "idx"))
        .def("NumRings", &RingInfo::numRings, python::args("self"))
        .def("IsRingFused", &RingInfo::isRingFused,
             python::args("self", "ringIdx"))
        .def("AreRingsFused", &RingInfo::areRingsFused,
             python::args("self", "ring1Idx", "ring2Idx"))
        .def("NumFusedBonds", &RingInfo::numFusedBonds,
             python::args("self", "ringIdx"))
        .def("AtomRings", atomRings, python::args("self"))
        .def("BondRings", bondRings, python::args("self"))
        .def("AtomMembers", atomMembers, python::args("self", "idx"))
        .def("BondMembers", bondMembers, python::args("self", "idx"))
        .def("AtomRingSizes", atomRingSizes, python::args("self", "idx"))
        .def("BondRingSizes", bondRingSizes, python::args("self", "idx"))
        .def("NumRingFamilies", &RingInfo::numRingFamilies,
             python::args("self"))
        .def("NumRelevantCycles", &RingInfo::numRelevantCycles,
             python::args("self"))
        .def("AtomRingFamilies", atomRingFamilies, python::args("self"))
        .def("BondRingFamilies", bondRingFamilies, python::args("self"))
        .def("AreRingFamiliesInitialized",
             &RingInfo::areRingFamiliesInitialized, python::args("self"))
        .def("AddRing", addRing,
             (python::arg("self"), python::arg("atomIds"),
              python::arg("bondIds")),
             "Adds a ring to the set. Be very careful with this operation.");
  };
};
}  // namespace RDKit
void wrap_ringinfo() { RDKit::ringinfo_wrapper::wrap(); }
