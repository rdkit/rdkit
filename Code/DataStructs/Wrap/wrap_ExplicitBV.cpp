//
//  Copyright (C) 2003-2021 greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDBoost/Wrap.h>
#include <DataStructs/BitVects.h>
#include <RDBoost/PySequenceHolder.h>
#include "wrap_helpers.h"

namespace python = boost::python;

// allows BitVects to be pickled
struct ebv_pickle_suite : rdkit_pickle_suite {
  static python::tuple getinitargs(const ExplicitBitVect &self) {
    std::string res = self.toString();
    python::object retval = python::object(
        python::handle<>(PyBytes_FromStringAndSize(res.c_str(), res.length())));
    return python::make_tuple(retval);
  };
};

python::list ExplicitToList(const ExplicitBitVect &sv) {
  python::list l;
  if (!sv.dp_bits) {
    return l;
  }

  auto count = sv.getNumBits();
  if (count) {
    l.append(0);
    l *= count;
    auto pos = sv.dp_bits->find_first();
    while (pos != boost::dynamic_bitset<>::npos) {
      l[pos] = 1;
      pos = sv.dp_bits->find_next(pos);
    }
  }
  return l;
}

std::string ebvClassDoc =
    "A class to store explicit bit vectors.\n\
\n\
This class is most useful for situations where the size of the vector\n\
is relatively small (tens of thousands or smaller).\n\
\n\
For larger vectors, use the _SparseBitVect_ class instead.\n\
\n\
As you would expect, _ExplicitBitVects_ support a set of binary operations\n\
so you can do things like:\n\
  bv3 = bv1 & bv2  (bitwise and)\n\
  bv3 = bv1 | bv2  (bitwise or)\n\
  bv3 = bv1 ^ bv2  (bitwise xor)\n\
  bv3 = ~bv1       (bitwise negation)\n\
\n\
Bits can be set and read using either the Set/UnsetBit() and GetBit() methods\n\
or by indexing (i.e. bv[i] = 1 or if bv[i]).\n\
\n";

struct EBV_wrapper {
  static void wrap() {
    python::class_<ExplicitBitVect, boost::shared_ptr<ExplicitBitVect>>(
        "ExplicitBitVect", ebvClassDoc.c_str(),
        python::init<unsigned int>(python::args("self", "size")))
        .def(python::init<std::string>(python::args("self", "pkl")))
        .def(python::init<unsigned int, bool>(
            python::args("self", "size", "bitsSet")))
        .def("SetBit", (bool(EBV::*)(unsigned int)) & EBV::setBit,
             python::args("self", "which"),
             "Turns on a particular bit.  Returns the original state of the "
             "bit.\n")
        .def("SetBitsFromList",
             (void (*)(EBV *, python::object))SetBitsFromList,
             python::args("self", "onBitList"),
             "Turns on a set of bits.  The argument should be a tuple or list "
             "of bit ids.\n")
        .def("UnSetBit", (bool(EBV::*)(unsigned int)) & EBV::unsetBit,
             python::args("self", "which"),
             "Turns off a particular bit.  Returns the original state of the "
             "bit.\n")
        .def("UnSetBitsFromList",
             (void (*)(EBV *, python::object))UnSetBitsFromList,
             python::args("self", "offBitList"),
             "Turns off a set of bits.  The argument should be a tuple or list "
             "of bit ids.\n")
        .def("GetBit", (bool(EBV::*)(unsigned int) const) & EBV::getBit,
             python::args("self", "which"), "Returns the value of a bit.\n")
        .def("GetNumBits", &EBV::getNumBits, python::args("self"),
             "Returns the number of bits in the vector (the vector's size).\n")
        .def("__len__", &EBV::getNumBits, python::args("self"))
        .def("GetNumOnBits", &EBV::getNumOnBits, python::args("self"),
             "Returns the number of on bits.\n")
        .def("GetNumOffBits", &EBV::getNumOffBits, python::args("self"),
             "Returns the number of off bits.\n")
        .def("__getitem__", (int (*)(const EBV &, int))get_VectItem,
             python::args("self", "which"))
        .def("__setitem__", (int (*)(EBV &, int, int))set_VectItem,
             python::args("self", "which", "val"))
        .def("GetOnBits", (IntVect(*)(const EBV &))GetOnBits,
             python::args("self"),
             "Returns a tuple containing IDs of the on bits.\n")
        .def("ToBinary", (python::object(*)(const EBV &))BVToBinary,
             python::args("self"),
             "Returns an internal binary representation of the vector.\n")
        .def("FromBase64", (void (*)(EBV &, const std::string &))InitFromBase64,
             python::args("self", "inD"),
             "Initializes the vector from a base64 encoded binary string.\n")
        .def("ToBase64", (std::string(*)(EBV &))ToBase64, python::args("self"),
             "Converts the vector to a base64 string (the base64 encoded "
             "version of the results of ToString()).\n")
        .def("ToList", (python::list(*)(const EBV &))ExplicitToList,
             python::args("self"),
             "Return the Bitvector as a python list (faster than list(vect))")
        .def(python::self & python::self)
        .def(python::self | python::self)
        .def(python::self ^ python::self)
        .def(python::self + python::self)
        .def(~python::self)
        .def(python::self == python::self)
        .def(python::self != python::self)
        .def(python::self += python::self)

        .def_pickle(ebv_pickle_suite());
  }
};

void wrap_EBV() { EBV_wrapper::wrap(); }
