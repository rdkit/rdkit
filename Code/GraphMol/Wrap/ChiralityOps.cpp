//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/python.h>

#include <string>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Chirality.h>

#include <RDBoost/Wrap.h>

namespace python = boost::python;
namespace RDKit {
struct chiralityops_wrapper {
  static void wrap() {
    RegisterVectorConverter<Chirality::StereoInfo>();

    python::def("FindPotentialStereo", &Chirality::findPotentialStereo,
                python::with_custodian_and_ward_postcall<0, 1>(),
                "find potential stereo elements in a molecule");
  };
};
}  // namespace RDKit

void wrap_chiralityops() { RDKit::chiralityops_wrapper::wrap(); }
