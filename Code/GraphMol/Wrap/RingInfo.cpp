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
  VECT_INT_VECT rings = self->atomRings();
  for (auto &ring : rings) {
    res.append(python::tuple(ring));
  }
  return python::tuple(res);
}
python::object bondRings(const RingInfo *self) {
  python::list res;
  VECT_INT_VECT rings = self->bondRings();
  for (auto &ring : rings) {
    res.append(python::tuple(ring));
  }
  return python::tuple(res);
}
void addRing(RingInfo *self,python::object atomRing, python::object bondRing){
  unsigned int nAts = python::extract<unsigned int>(atomRing.attr("__len__")());
  unsigned int nBnds = python::extract<unsigned int>(bondRing.attr("__len__")());
  if(nAts != nBnds) throw_value_error("list sizes must match");
  if(!self->isInitialized()) self->initialize();
  INT_VECT aring(nAts);
  INT_VECT bring(nAts);
  for (unsigned int i = 0; i < nAts; ++i) {
    aring[i] = python::extract<int>(atomRing[i])();
    bring[i] = python::extract<int>(bondRing[i])();
  }
  self->addRing(aring,bring);
}
}

namespace RDKit {
std::string classDoc = "contains information about a molecule's rings\n";

struct ringinfo_wrapper {
  static void wrap() {
    python::class_<RingInfo>("RingInfo", classDoc.c_str(), python::no_init)
        .def("IsAtomInRingOfSize", &RingInfo::isAtomInRingOfSize)
        .def("MinAtomRingSize", &RingInfo::minAtomRingSize)
        .def("IsBondInRingOfSize", &RingInfo::isBondInRingOfSize)
        .def("MinBondRingSize", &RingInfo::minBondRingSize)
        .def("NumAtomRings", &RingInfo::numAtomRings)
        .def("NumBondRings", &RingInfo::numBondRings)
        .def("NumRings", &RingInfo::numRings)
        .def("AtomRings", atomRings)
        .def("BondRings", bondRings)
        .def("AddRing", addRing, (python::arg("self"),python::arg("atomIds"),python::arg("bondIds")),
         "Adds a ring to the set. Be very careful with this operation.");
  };
};
}  // end of namespace
void wrap_ringinfo() { RDKit::ringinfo_wrapper::wrap(); }
