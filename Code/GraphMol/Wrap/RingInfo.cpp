// $Id$
//
//  Created by Greg Landrum: January 2007
//
//   @@ All Rights Reserved  @@
//

#define NO_IMPORT_ARRAY
#include <boost/python.hpp>

#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>

namespace python = boost::python;


namespace {
  using namespace RDKit;
  python::object atomRings(const RingInfo *self){
    python::list res;
    VECT_INT_VECT rings=self->atomRings();
    for(VECT_INT_VECT_I ringIt=rings.begin();ringIt!=rings.end();++ringIt){
      res.append(python::tuple(*ringIt));
    }
    return python::tuple(res);
  }
  python::object bondRings(const RingInfo *self){
    python::list res;
    VECT_INT_VECT rings=self->bondRings();
    for(VECT_INT_VECT_I ringIt=rings.begin();ringIt!=rings.end();++ringIt){
      res.append(python::tuple(*ringIt));
    }
    return python::tuple(res);
  }
}

namespace RDKit{
    std::string classDoc="contains information about a molecule's rings\n";
  
struct ringinfo_wrapper {
  static void wrap(){
    python::class_<RingInfo>("RingInfo",classDoc.c_str(),python::no_init)
      .def("IsAtomInRingOfSize",&RingInfo::isAtomInRingOfSize)
      .def("IsBondInRingOfSize",&RingInfo::isBondInRingOfSize)
      .def("IsBondInRingOfSize",&RingInfo::isBondInRingOfSize)
      .def("NumAtomRings",&RingInfo::numAtomRings)
      .def("NumBondRings",&RingInfo::numBondRings)
      .def("NumRings",&RingInfo::numRings)
      .def("AtomRings",atomRings)
      .def("BondRings",bondRings)
      ;
  };
};
}// end of namespace
void wrap_ringinfo() {
  RDKit::ringinfo_wrapper::wrap();
}
