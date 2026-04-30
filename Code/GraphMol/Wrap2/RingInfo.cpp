//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>

#include <RDBoost/Wrap_nb.h>

#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace {
using namespace RDKit;
nb::tuple atomRings(const RingInfo *self) {
  nb::list res;
  for (const auto &ring : self->atomRings()) {
    nb::list ringAsList;
    for (const auto idx : ring) {
      ringAsList.append(idx);
    }
    res.append(nb::tuple(ringAsList));
  }
  return nb::tuple(res);
}
nb::tuple bondRings(const RingInfo *self) {
  nb::list res;
  for (const auto &ring : self->bondRings()) {
    nb::list ringAsList;
    for (const auto idx : ring) {
      ringAsList.append(idx);
    }
    res.append(nb::tuple(ringAsList));
  }
  return nb::tuple(res);
}
nb::tuple atomMembers(const RingInfo *self, unsigned int idx) {
  nb::list res;
  for (const auto v : self->atomMembers(idx)) {
    res.append(v);
  }
  return nb::tuple(res);
}
nb::tuple bondMembers(const RingInfo *self, unsigned int idx) {
  nb::list res;
  for (const auto v : self->bondMembers(idx)) {
    res.append(v);
  }
  return nb::tuple(res);
}
nb::tuple atomRingSizes(const RingInfo *self, unsigned int idx) {
  nb::list res;
  for (const auto v : self->atomRingSizes(idx)) {
    res.append(v);
  }
  return nb::tuple(res);
}
nb::tuple bondRingSizes(const RingInfo *self, unsigned int idx) {
  nb::list res;
  for (const auto v : self->bondRingSizes(idx)) {
    res.append(v);
  }
  return nb::tuple(res);
}

nb::tuple atomRingFamilies(const RingInfo *self) {
  nb::list res;
  for (const auto &ring : self->atomRingFamilies()) {
    nb::list ringAsList;
    for (const auto idx : ring) {
      ringAsList.append(idx);
    }
    res.append(nb::tuple(ringAsList));
  }
  return nb::tuple(res);
}
nb::tuple bondRingFamilies(const RingInfo *self) {
  nb::list res;
  for (const auto &ring : self->bondRingFamilies()) {
    nb::list ringAsList;
    for (const auto idx : ring) {
      ringAsList.append(idx);
    }
    res.append(nb::tuple(ringAsList));
  }
  return nb::tuple(res);
}

void addRing(RingInfo *self, const nb::object &atomRing,
             const nb::object &bondRing) {
  auto atomIds = pythonObjectToVect<int>(atomRing);
  auto bondIds = pythonObjectToVect<int>(bondRing);
  const unsigned int nAts = atomIds ? atomIds->size() : 0;
  const unsigned int nBnds = bondIds ? bondIds->size() : 0;
  if (nAts != nBnds) {
    throw ValueErrorException("list sizes must match");
  }
  if (!self->isInitialized()) {
    self->initialize();
  }
  INT_VECT aring(atomIds->begin(), atomIds->end());
  INT_VECT bring(bondIds->begin(), bondIds->end());
  self->addRing(aring, bring);
}
}  // namespace

namespace RDKit {
std::string classDoc =
    R"DOC(contains information about a molecule's rings
)DOC";

struct ringinfo_wrapper {
  static void wrap(nb::module_ &m) {
    nb::class_<RingInfo>(m, "RingInfo", classDoc.c_str())
        .def("IsAtomInRingOfSize", &RingInfo::isAtomInRingOfSize, "idx"_a,
             "size"_a)
        .def("MinAtomRingSize", &RingInfo::minAtomRingSize, "idx"_a)
        .def("AreAtomsInSameRing", &RingInfo::areAtomsInSameRing, "idx1"_a,
             "idx2"_a)
        .def("AreAtomsInSameRingOfSize", &RingInfo::areAtomsInSameRingOfSize,
             "idx1"_a, "idx2"_a, "size"_a)
        .def("IsBondInRingOfSize", &RingInfo::isBondInRingOfSize, "idx"_a,
             "size"_a)
        .def("MinBondRingSize", &RingInfo::minBondRingSize, "idx"_a)
        .def("AreBondsInSameRing", &RingInfo::areBondsInSameRing, "idx1"_a,
             "idx2"_a)
        .def("AreBondsInSameRingOfSize", &RingInfo::areBondsInSameRingOfSize,
             "idx1"_a, "idx2"_a, "size"_a)
        .def("NumAtomRings", &RingInfo::numAtomRings, "idx"_a)
        .def("NumBondRings", &RingInfo::numBondRings, "idx"_a)
        .def("NumRings", &RingInfo::numRings)
        .def("IsRingFused", &RingInfo::isRingFused, "ringIdx"_a)
        .def("AreRingsFused", &RingInfo::areRingsFused, "ring1Idx"_a,
             "ring2Idx"_a)
        .def("NumFusedBonds", &RingInfo::numFusedBonds, "ringIdx"_a)
        .def("AtomRings", atomRings)
        .def("BondRings", bondRings)
        .def("AtomMembers", atomMembers, "idx"_a)
        .def("BondMembers", bondMembers, "idx"_a)
        .def("AtomRingSizes", atomRingSizes, "idx"_a)
        .def("BondRingSizes", bondRingSizes, "idx"_a)
        .def("NumRingFamilies", &RingInfo::numRingFamilies)
        .def("NumRelevantCycles", &RingInfo::numRelevantCycles)
        .def("AtomRingFamilies", atomRingFamilies)
        .def("BondRingFamilies", bondRingFamilies)
        .def("AreRingFamiliesInitialized",
             &RingInfo::areRingFamiliesInitialized)
        .def(
            "AddRing", addRing, "atomIds"_a, "bondIds"_a,
            R"DOC(Adds a ring to the set. Be very careful with this operation.)DOC");
  };
};
}  // namespace RDKit
void wrap_ringinfo(nb::module_ &m) { RDKit::ringinfo_wrapper::wrap(m); }
